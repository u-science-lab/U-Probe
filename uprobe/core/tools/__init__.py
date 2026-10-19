from pathlib import Path
from .aligner import (
    build_bowtie2_index, build_bowtie_index, build_blast_db, build_mmseqs_index
)
from uprobe.core.gen.fun import extract_trans_seqs
from uprobe.core.utils import get_logger

log = get_logger(__name__)

_BOWTIE2_SMALL_SUFFIXES = (
    ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"
)
_BOWTIE2_LARGE_SUFFIXES = tuple(suffix + "l" for suffix in _BOWTIE2_SMALL_SUFFIXES)
_BOWTIE_SMALL_SUFFIXES = (
    ".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt"
)
_BOWTIE_LARGE_SUFFIXES = tuple(suffix + "l" for suffix in _BOWTIE_SMALL_SUFFIXES)
_BLAST_SUFFIXES = (".nin", ".nhr", ".nsq")
# mmseqs createdb writes <prefix>, <prefix>_h and their .dbtype/.index; createindex adds <prefix>.idx*.
_MMSEQS_SUFFIXES = ("", ".dbtype", ".index", "_h", "_h.dbtype", "_h.index", ".idx", ".idx.dbtype", ".idx.index")


def _blast_index_exists(index_prefix: Path) -> bool:
    return all(Path(str(index_prefix) + suffix).is_file()
               and Path(str(index_prefix) + suffix).stat().st_size > 0
               for suffix in _BLAST_SUFFIXES)


def _bowtie2_index_exists(index_prefix: Path) -> bool:
    """Return True for either a complete .bt2 or .bt2l index."""
    prefix = str(index_prefix)
    return (
        all(Path(prefix + suffix).is_file() and Path(prefix + suffix).stat().st_size > 0
            for suffix in _BOWTIE2_SMALL_SUFFIXES)
        or all(Path(prefix + suffix).is_file() and Path(prefix + suffix).stat().st_size > 0
               for suffix in _BOWTIE2_LARGE_SUFFIXES)
    )


def _bowtie_index_exists(index_prefix: Path) -> bool:
    """Return True for either a complete .ebwt or .ebwtl (Bowtie 1) index."""
    prefix = str(index_prefix)
    return any(
        all(Path(prefix + suffix).is_file() and Path(prefix + suffix).stat().st_size > 0
            for suffix in suffixes)
        for suffixes in (_BOWTIE_SMALL_SUFFIXES, _BOWTIE_LARGE_SUFFIXES)
    )


def _mmseqs_index_exists(index_prefix: Path) -> bool:
    return all(Path(str(index_prefix) + suffix).is_file()
               and Path(str(index_prefix) + suffix).stat().st_size > 0
               for suffix in _MMSEQS_SUFFIXES)


_GENOME_INDEX_EXISTS = {
    "bowtie2": _bowtie2_index_exists,
    "blast": _blast_index_exists,
    "mmseqs": _mmseqs_index_exists,
}
_TRANSCRIPT_INDEX_EXISTS = {
    "bowtie2": _bowtie2_index_exists,
    "bowtie": _bowtie_index_exists,
    "blast": _blast_index_exists,
}
# Bowtie 1 only serves mapped_genes (aligner: bowtie), which searches the
# transcriptome, so it gets no genome index.
_TRANSCRIPT_ONLY = ("bowtie",)


def transcript_fasta(fasta: Path, outdir: Path) -> Path:
    """The transcript.fa a transcriptome index in ``outdir`` is built from.

    It lies next to the genome FASTA or in ``outdir``; when neither exists,
    the path next to the genome FASTA, where it is extracted.
    """
    candidates = (fasta.parent / "transcript.fa", outdir / "transcript.fa")
    return next((p for p in candidates if p.is_file() and p.stat().st_size > 0), candidates[0])


def build_transcripts_index(gtf: Path,
                             fasta: Path,
                             outdir: Path,
                             threads: int = 10,
                             aligner: str = "bowtie2"
                            ) -> str:
    if aligner not in _TRANSCRIPT_INDEX_EXISTS:
        raise ValueError(f"Unsupported transcript aligner: {aligner}")
    expected_dir = f"{aligner}_transcript"
    if outdir.name in {f"{name}_transcript" for name in _TRANSCRIPT_INDEX_EXISTS} and outdir.name != expected_dir:
        raise ValueError(f"{aligner} transcript index must use {expected_dir}, got {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    index_prefix = outdir / fasta.stem
    index_exists = _TRANSCRIPT_INDEX_EXISTS[aligner]
    if index_exists(index_prefix):
        return str(index_prefix)
    trans_fasta_path = transcript_fasta(fasta, outdir)
    if not trans_fasta_path.is_file() or trans_fasta_path.stat().st_size == 0:
        extract_trans_seqs(gtf, fasta, trans_fasta_path)
    if not trans_fasta_path.is_file() or trans_fasta_path.stat().st_size == 0:
        raise ValueError(f"Transcript FASTA is empty or missing: {trans_fasta_path}")
    if aligner == "bowtie2":
        build_bowtie2_index(trans_fasta_path, index_prefix, threads)
    elif aligner == "bowtie":
        build_bowtie_index(trans_fasta_path, index_prefix, threads)
    else:
        build_blast_db(trans_fasta_path, index_prefix, title=f"{fasta.stem}_transcript")
    if not index_exists(index_prefix):
        raise RuntimeError(f"Incomplete {aligner} transcript index: {index_prefix}")
    return str(index_prefix)

def build_genome(genome: dict,
                 threads: int = 10
                 ) -> dict:
    """Build the genome index using the provided fasta file.

    An index is skipped only when all of its files exist, so a build that was
    interrupted (leaving the <aligner>_genome folder behind) is redone.
    """
    fasta_path = Path(genome['fasta'])
    prefix = fasta_path.stem
    aligner_index: list = genome['align_index']
    for aligner in aligner_index:
        if aligner in _TRANSCRIPT_ONLY:
            continue
        if aligner not in _GENOME_INDEX_EXISTS:
            raise NotImplementedError(f"aligner {aligner} is not implemented")
        index_dir = fasta_path.parent / f"{aligner}_genome"/ prefix
        if _GENOME_INDEX_EXISTS[aligner](index_dir):
            log.info(f"index already exists: {index_dir}")
            continue
        log.info(f"index does not exist or is incomplete, building {aligner} index for {prefix}: {index_dir}")
        index_dir.parent.mkdir(parents=True, exist_ok=True)
        if aligner == "bowtie2":
            build_bowtie2_index(fasta_path, index_dir, threads)
        elif aligner == "blast":
            build_blast_db(fasta_path, index_dir, title=prefix)
        else:
            build_mmseqs_index(fasta_path, str(index_dir))
    for aligner in aligner_index:
        if aligner not in _TRANSCRIPT_INDEX_EXISTS:
            continue
        tran_index_dir = fasta_path.parent / f"{aligner}_transcript"
        build_transcripts_index(gtf=Path(genome['gtf']),
                                fasta=fasta_path,
                                outdir=tran_index_dir,
                                threads=threads,
                                aligner=aligner)
    return genome
