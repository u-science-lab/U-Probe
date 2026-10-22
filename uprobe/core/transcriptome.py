"""Register a transcriptome FASTA and a transcript table as a genomes.yaml entry.

``uprobe genome add --transcriptome transcripts.fa --records transcripts.tsv`` for designs
that start from transcript sequences instead of a genome and its GTF. Each FASTA record is a
transcript, written as a single-exon, plus-strand "contig" of a synthetic GTF whose gene names,
gene ids and tags (MANE_Select, Ensembl_canonical) come from the table, so ``source: mRNA`` with
``transcripts: mane`` or ``canonical`` and the transcriptome off-target search work as with a
genome. The entry keeps everything it derives (the GTF, transcript.fa with the
``<gene_id>_<transcript_id>`` record names the off-target search reads, and the indexes) in its own
directory (``index_dir``), not next to the FASTA, and can point at indexes that already exist.
"""
import os
import re
import typing as T
from pathlib import Path

from .genome_fetch import _NAME_RE, default_align_index, load_genomes, update_genomes_yaml
from .genome_paths import KINDS
from .utils import get_logger

log = get_logger(__name__)

TRANSCRIPT_COLUMNS = ('transcript_id', 'transcript', 'ensembl_transcript_id', 'enst', 'refseq_transcript_id', 'id')
GENE_COLUMNS = ('gene_name', 'gene', 'gene_symbol', 'symbol', 'hgnc_symbol', 'external_gene_name')
GENE_ID_COLUMNS = ('gene_id', 'ensembl_gene_id', 'ensg')
TAG_COLUMNS = ('tag', 'tags')
MANE_COLUMNS = ('mane_select', 'mane', 'is_mane_select', 'mane_status')
CANONICAL_COLUMNS = ('ensembl_canonical', 'canonical', 'is_canonical', 'is_ensembl_canonical')
CDS_COLUMNS = (('cds_start', 'cds_end'), ('cds_start_in_transcript', 'cds_end_in_transcript'))
TRUE_WORDS = {'1', 'true', 't', 'yes', 'y'}
MANE_TAG, CANONICAL_TAG = 'MANE_Select', 'Ensembl_canonical'
_VERSION = re.compile(r'\.\d+$')


class TranscriptomeError(Exception):
    pass


def _find_column(columns: T.Sequence[str], aliases: T.Sequence[str], given: T.Optional[str], what: str,
                 required: bool = False) -> T.Optional[str]:
    lower = {str(c).strip().lower(): c for c in columns}
    if given is not None:
        if given not in columns:
            raise TranscriptomeError(f'the {what} column {given!r} is not in the table; its columns are: '
                                     f'{", ".join(map(str, columns))}')
        return given
    for alias in aliases:
        if alias in lower:
            return lower[alias]
    if required:
        raise TranscriptomeError(f'no {what} column in the table (looked for {", ".join(aliases)}); its columns '
                                 f'are: {", ".join(map(str, columns))}. Name it with the matching option.')
    return None


def _truthy(value, mane: bool = False) -> bool:
    text = str(value).strip().lower()
    if text in TRUE_WORDS:
        return True
    return 'select' in text if mane else text in ('canonical', 'ensembl_canonical')


def read_records(path, transcript_column=None, gene_column=None, gene_id_column=None, tag_column=None,
                 mane_column=None, canonical_column=None) -> T.Dict[str, dict]:
    """{transcript id: {gene_name, gene_id, tags, cds}} of a transcript table (TSV, or CSV when named .csv)."""
    from .gene_expression import read_table
    try:
        table = read_table(path)
    except Exception as exc:
        raise TranscriptomeError(f'cannot read the records table {path}: {exc}') from exc
    columns = list(table.columns)
    transcript = _find_column(columns, TRANSCRIPT_COLUMNS, transcript_column, 'transcript id', required=True)
    gene = _find_column(columns, GENE_COLUMNS, gene_column, 'gene name')
    gene_id = _find_column(columns, GENE_ID_COLUMNS, gene_id_column, 'gene id')
    if gene is None and gene_id is None:
        raise TranscriptomeError('no gene column in the table (looked for '
                                 f'{", ".join(GENE_COLUMNS + GENE_ID_COLUMNS)}); its columns are: '
                                 f'{", ".join(map(str, columns))}')
    tag = _find_column(columns, TAG_COLUMNS, tag_column, 'tag')
    mane = _find_column(columns, MANE_COLUMNS, mane_column, 'MANE Select')
    canonical = _find_column(columns, CANONICAL_COLUMNS, canonical_column, 'Ensembl canonical')
    cds = next(((a, b) for a, b in CDS_COLUMNS
                if a in {str(c).lower() for c in columns} and b in {str(c).lower() for c in columns}), None)
    if cds:
        by_lower = {str(c).lower(): c for c in columns}
        cds = (by_lower[cds[0]], by_lower[cds[1]])
    records: T.Dict[str, dict] = {}
    for row in table.to_dict('records'):
        name = row[transcript]
        if not isinstance(name, str) or not name.strip():
            continue
        name = name.strip()
        gene_name = row[gene] if gene and isinstance(row[gene], str) and row[gene].strip() else None
        gid = row[gene_id] if gene_id and isinstance(row[gene_id], str) and row[gene_id].strip() else None
        tags = set()
        if tag and isinstance(row[tag], str):
            tags.update(t.strip() for t in re.split(r'[,;|]', row[tag]) if t.strip())
        if mane and isinstance(row[mane], str) and _truthy(row[mane], mane=True):
            tags.add(MANE_TAG)
        if canonical and isinstance(row[canonical], str) and _truthy(row[canonical]):
            tags.add(CANONICAL_TAG)
        span = None
        if cds:
            try:
                span = (int(float(row[cds[0]])), int(float(row[cds[1]])))
            except (TypeError, ValueError):
                span = None
        entry = records.setdefault(name, {'gene_name': gene_name or gid, 'gene_id': gid or gene_name,
                                          'tags': set(), 'cds': span})
        entry['tags'] |= tags
    if not records:
        raise TranscriptomeError(f'the records table {path} lists no transcript')
    return records


def _fasta_records(path: Path) -> T.Iterator[T.Tuple[str, T.Iterator[str]]]:
    """(name, sequence lines) per record; the name is the header up to the first whitespace."""
    with open(path, 'r') as handle:
        name, lines = None, []
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    yield name, iter(lines)
                name, lines = line[1:].split()[0] if line[1:].strip() else '', []
            elif name is not None:
                lines.append(line.rstrip('\r\n'))
        if name is not None:
            yield name, iter(lines)


def _safe(identifier: str) -> str:
    """A gene id as the off-target search reads it: record names are <gene_id>_<transcript_id>, cut at the first '_'."""
    return identifier.replace('_', '-')


def _gtf_lines(seqid, length, record, transcript_id, gene_id) -> T.List[str]:
    attrs = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; gene_name "{record["gene_name"]}";'
    attrs += ''.join(f' tag "{tag}";' for tag in sorted(record['tags']))
    lines = [f'{seqid}\tuprobe\ttranscript\t1\t{length}\t.\t+\t.\t{attrs}',
             f'{seqid}\tuprobe\texon\t1\t{length}\t.\t+\t.\t{attrs}']
    if record.get('cds'):
        start, end = record['cds']
        if 1 <= start <= end <= length:
            lines.append(f'{seqid}\tuprobe\tCDS\t{start}\t{end}\t.\t+\t0\t{attrs}')
    return lines


def add_transcriptome(name: str, fasta: T.Union[str, Path], records: T.Union[str, Path],
                      out_dir: T.Union[str, Path] = 'genomes', genomes_yaml: T.Union[str, Path] = 'genomes.yaml',
                      *, index_dir: T.Optional[T.Union[str, Path]] = None,
                      existing_indexes: T.Optional[T.Dict[str, str]] = None,
                      align_index: T.Optional[T.List[str]] = None, jellyfish: T.Optional[bool] = None,
                      unlisted: str = 'skip', species: T.Optional[str] = None, dry_run: bool = False,
                      **columns) -> dict:
    """Write ``out_dir/<name>/annotation.gtf`` and ``transcript.fa`` and register the entry.

    ``columns`` name table columns when the defaults are not found: transcript_column, gene_column,
    gene_id_column, tag_column, mane_column, canonical_column. ``unlisted`` is what to do with FASTA
    records that the table lacks: 'skip' (left out of the annotation and of the off-target index) or
    'include' (each is its own gene, named by its record). ``existing_indexes`` maps <aligner>_<kind>
    (e.g. bowtie_transcript) to the prefix of an index that exists. Returns {name, entry, summary,
    directory, genomes_yaml, changed}.
    """
    if not _NAME_RE.fullmatch(name):
        raise TranscriptomeError(f"Invalid genome name '{name}': use letters, digits, '.', '_' or '-'.")
    if unlisted not in ('skip', 'include'):
        raise TranscriptomeError(f"unlisted must be skip or include, got {unlisted!r}")
    fasta, records = Path(fasta).resolve(), Path(records)
    if not fasta.is_file():
        raise TranscriptomeError(f'transcriptome FASTA not found: {fasta}')
    if fasta.suffix == '.gz':
        raise TranscriptomeError(f'{fasta.name} is compressed; decompress it first (the aligners and pyfaidx need '
                                 'a plain FASTA)')
    table = read_records(records, **{key: value for key, value in columns.items() if value is not None})
    by_version = {}
    for key in table:
        by_version.setdefault(_VERSION.sub('', key), []).append(key)
    for kind_name in (existing_indexes or {}):
        aligner, _, kind = kind_name.partition('_')
        if kind not in KINDS or not aligner:
            raise TranscriptomeError(f"existing index {kind_name!r}: use <aligner>_<kind> with kind genome or "
                                     "transcript, e.g. bowtie_transcript")

    directory = (Path(out_dir) / name).resolve()
    root = Path(index_dir).resolve() if index_dir else directory
    genomes_yaml = Path(genomes_yaml).resolve()
    gtf_path, transcript_path, link = directory / 'annotation.gtf', root / 'transcript.fa', directory / 'transcriptome.fa'
    summary = {'fasta_records': 0, 'transcripts': 0, 'genes': 0, 'unlisted': 0, 'table_only': 0,
               'gene_ids_changed': 0, 'with_cds': 0}
    if not dry_run:
        directory.mkdir(parents=True, exist_ok=True)
        root.mkdir(parents=True, exist_ok=True)
    genes, seen, matched = set(), set(), set()
    out = None if dry_run else open(str(transcript_path) + '.tmp', 'w')
    gtf_out = None if dry_run else open(str(gtf_path) + '.tmp', 'w')
    try:
        for record_name, lines in _fasta_records(fasta):
            summary['fasta_records'] += 1
            if record_name in seen:
                raise TranscriptomeError(f'the FASTA has two records named {record_name!r}')
            seen.add(record_name)
            key = record_name if record_name in table else next(
                (k for k in by_version.get(_VERSION.sub('', record_name), []) if len(by_version[_VERSION.sub('', record_name)]) == 1), None)
            if key is None:
                summary['unlisted'] += 1
                if unlisted == 'skip':
                    continue
                entry = {'gene_name': record_name, 'gene_id': record_name, 'tags': set(), 'cds': None}
            else:
                entry = table[key]
                matched.add(key)
            sequence = ''.join(lines)
            gene_id = _safe(entry['gene_id'])
            summary['gene_ids_changed'] += gene_id != entry['gene_id']
            summary['with_cds'] += bool(entry.get('cds')) and 1 <= entry['cds'][0] <= entry['cds'][1] <= len(sequence)
            if gtf_out:
                gtf_out.write('\n'.join(_gtf_lines(record_name, len(sequence), entry, record_name, gene_id)) + '\n')
            genes.add(gene_id)
            summary['transcripts'] += 1
            if out:
                out.write(f'>{gene_id}_{record_name}\n{sequence}\n')
    finally:
        for handle in (out, gtf_out):
            if handle:
                handle.close()
    summary['table_only'] = len(set(table) - matched)
    summary['genes'] = len(genes)
    if not summary['transcripts']:
        if out:
            os.unlink(str(transcript_path) + '.tmp')
            os.unlink(str(gtf_path) + '.tmp')
        raise TranscriptomeError(f'none of the {summary["fasta_records"]} FASTA records is in the records table; '
                                 'check that the transcript column holds the FASTA record names (the header up to '
                                 'the first space)')
    if summary['unlisted'] and unlisted == 'skip':
        log.warning(f"{summary['unlisted']} of {summary['fasta_records']} FASTA records are not in the records "
                    "table and are left out of the annotation and the off-target index (--unlisted include keeps "
                    "them as genes of their own)")
    if summary['table_only']:
        log.warning(f"{summary['table_only']} transcripts of the records table are not in the FASTA")
    if summary['gene_ids_changed']:
        log.warning(f"{summary['gene_ids_changed']} gene ids contain '_', which the off-target search reads as the "
                    "end of the gene id; written with '-' instead")

    entry = {
        'description': f'Transcriptome {fasta.name} ({summary["transcripts"]:,} transcripts, {summary["genes"]:,} '
                       'genes; single-exon synthetic annotation)',
        'fasta': str(link), 'gtf': str(gtf_path), 'out': str(directory),
        'align_index': list(align_index) if align_index is not None else default_align_index(),
        'jellyfish': bool(jellyfish) if jellyfish is not None else False,
        'source': {'provider': 'uprobe genome add', 'transcriptome': str(fasta), 'records': str(Path(records).resolve()),
                   'transcripts': summary['transcripts'], 'genes': summary['genes']},
    }
    if species:
        entry['species'] = species
    if index_dir:
        entry['index_dir'] = str(root)
    else:
        entry['index_dir'] = str(directory)
    if existing_indexes:
        entry['indexes'] = {key: str(Path(value).resolve()) for key, value in existing_indexes.items()}
    result = {'name': name, 'entry': entry, 'summary': summary, 'directory': directory,
              'genomes_yaml': genomes_yaml, 'changed': False}
    if dry_run:
        return result
    os.replace(str(gtf_path) + '.tmp', gtf_path)
    os.replace(str(transcript_path) + '.tmp', transcript_path)
    if link.is_symlink() or link.exists():
        link.unlink()
    try:
        link.symlink_to(fasta)
    except OSError:
        import shutil
        shutil.copyfile(fasta, link)
    for key, value in (existing_indexes or {}).items():
        from .tools import _GENOME_INDEX_EXISTS, _TRANSCRIPT_INDEX_EXISTS
        aligner, _, kind = key.partition('_')
        exists = (_GENOME_INDEX_EXISTS if kind == 'genome' else _TRANSCRIPT_INDEX_EXISTS).get(aligner)
        if exists is None or not exists(Path(value)):
            log.warning(f'indexes.{key}: no complete {aligner} index found at the prefix {value}; runs that '
                        'search it will fail until it exists')
    result['changed'] = update_genomes_yaml(genomes_yaml, name, entry)
    log.info(f"Transcriptome '{name}' is registered in {genomes_yaml}")
    return result
