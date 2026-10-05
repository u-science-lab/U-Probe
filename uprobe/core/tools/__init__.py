from pathlib import Path
from .aligner import (
    build_bowtie2_index, build_blast_db, build_mmseqs_index
)
from uprobe.core.gen.fun import extract_trans_seqs
from uprobe.core.utils import get_logger

log = get_logger(__name__)

_BOWTIE2_SMALL_SUFFIXES = (
    ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"
)
_BOWTIE2_LARGE_SUFFIXES = tuple(suffix + "l" for suffix in _BOWTIE2_SMALL_SUFFIXES)
_BLAST_SUFFIXES = (".nin", ".nhr", ".nsq")


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

def build_transcripts_index(gtf: Path,
                             fasta: Path,
                             outdir: Path,
                             threads: int = 10,
                             aligner: str = "bowtie2"
                            ) -> str:
    if aligner not in {"bowtie2", "blast"}:
        raise ValueError(f"Unsupported transcript aligner: {aligner}")
    expected_dir = f"{aligner}_transcript"
    if outdir.name in {"blast_transcript", "bowtie2_transcript"} and outdir.name != expected_dir:
        raise ValueError(f"{aligner} transcript index must use {expected_dir}, got {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    index_prefix = outdir / fasta.stem
    index_exists = _bowtie2_index_exists if aligner == "bowtie2" else _blast_index_exists
    if index_exists(index_prefix):
        return str(index_prefix)
    candidates = (fasta.parent / "transcript.fa", outdir / "transcript.fa")
    trans_fasta_path = next((p for p in candidates if p.is_file() and p.stat().st_size > 0),
                           candidates[0])
    if not trans_fasta_path.is_file() or trans_fasta_path.stat().st_size == 0:
        extract_trans_seqs(gtf, fasta, trans_fasta_path)
    if not trans_fasta_path.is_file() or trans_fasta_path.stat().st_size == 0:
        raise ValueError(f"Transcript FASTA is empty or missing: {trans_fasta_path}")
    if aligner == "bowtie2":
        build_bowtie2_index(trans_fasta_path, index_prefix, threads)
    else:
        build_blast_db(trans_fasta_path, index_prefix, title=f"{fasta.stem}_transcript")
    if not index_exists(index_prefix):
        raise RuntimeError(f"Incomplete {aligner} transcript index: {index_prefix}")
    return str(index_prefix)

def build_genome(genome: dict,
                 threads: int = 10
                 ) -> dict:
    """Build the genome index using the provided fasta file."""
    fasta_path = Path(genome['fasta'])
    prefix = fasta_path.stem

    bowtie2_index_files = [
                    f"{prefix}.1.bt2",
                    f"{prefix}.2.bt2",
                    f"{prefix}.3.bt2",
                    f"{prefix}.4.bt2",
                    f"{prefix}.rev.1.bt2",
                    f"{prefix}.rev.2.bt2"
                    ]
    blast_index_files = [
                    f"{prefix}.ndb",
                    f"{prefix}.nin",
                    f"{prefix}.nhr",
                    f"{prefix}.nsq"
                    ]
    mmseqs_index_files = [
                    f"{prefix}.db",
                    f"{prefix}.dbtype",
                    f"{prefix}.dbstat",
                    f"{prefix}.dbmeta"
                    ]
    aligner_index: list = genome['align_index']
    for aligner in aligner_index:
        log.info(f"building {aligner} index for {prefix}")
        index_dir = fasta_path.parent / f"{aligner}_genome"/ prefix
        if index_dir.parent.exists():
            log.info(f"index already exists: {index_dir.parent}")
        else:
            log.info(f"index does not exist, building it now: {index_dir.parent}")
            index_dir.parent.mkdir(parents=True, exist_ok=True)
            if aligner == "bowtie2":
                if all((index_dir.parent / file_name).exists() for file_name in bowtie2_index_files):
                    log.info(f"index already exists: {index_dir.parent}")
                else:
                    build_bowtie2_index(fasta_path, index_dir)
            elif aligner == "blast":
                if all((index_dir.parent / file_name).exists() for file_name in blast_index_files):
                    log.info(f"index already exists: {index_dir.parent}")
                else:
                    build_blast_db(fasta_path, index_dir, title=prefix)
            elif aligner == "mmseqs":
                if all((index_dir.parent / file_name).exists() for file_name in mmseqs_index_files):
                    log.info(f"index already exists: {index_dir.parent}")
                else:
                    build_mmseqs_index(fasta_path, str(index_dir))
            else:
                raise NotImplementedError(f"aligner {aligner} is not implemented")
    for aligner in aligner_index:
        if aligner not in {"bowtie2", "blast"}:
            continue
        tran_index_dir = fasta_path.parent / f"{aligner}_transcript"
        build_transcripts_index(gtf=Path(genome['gtf']),
                                fasta=fasta_path,
                                outdir=tran_index_dir,
                                threads=threads,
                                aligner=aligner)
    return genome
