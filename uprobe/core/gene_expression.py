"""Expression-aware off-target counting.

A protocol may supply a gene expression table for the sample the probes are
used in::

    gene_expression:
      table: hek293_tpm.tsv    # TSV, CSV (name contains .csv) or GTEx GCT; may be gzipped
      level: gene              # default; transcript: the table lists transcripts (see below)
      gene_column: gene_id     # default: the first column; a list of columns is also accepted
      value_column: TPM        # default: the only numeric column
      min_value: 1             # genes with a lower value are not off-targets
      min_fraction_of_target: 0.05  # optional: nor are genes below 5% of the target's value
      absent: count            # default; ignore: genes missing from the table are not off-targets

``mapped_genes`` attributes then ignore off-target genes whose value is below
``min_value``: an RNA that is not there cannot take up a probe. With
``min_fraction_of_target``, an off-target gene must also reach that fraction
of the target gene's value: the threshold for a target's off-targets is
max(min_value, fraction x target value), so a low-level pseudogene of a highly
expressed target does not count, as its share of the spots would be small.
A target absent from the table gets min_value alone (the run warns). Genes are
looked up by Ensembl gene id (with or without version) and then by gene name
(symbol), so a table keyed by either works. Off-target genes absent from the
table still count (conservative) unless ``absent: ignore``; either way the
run reports how many there were. The target gene always counts. Relative
paths are resolved from the directory where uprobe runs, as for genomes.yaml.

With ``level: transcript`` the table lists transcripts (salmon's quant.sf,
kallisto's abundance.tsv, RSEM's isoforms.results), keyed by transcript id
(``transcript_column``, default the first column; with or without version).
An off-target gene's value for a sequence is then the sum of the values of
the gene's transcripts that the sequence hits, not the whole gene's: false
signal comes only from RNA molecules that carry the site, and a gene whose
expression lies on isoforms without it does not count. The genome annotation,
from which the searched transcriptome is built, assigns transcripts to genes;
a target's value is the sum over all its transcripts in the table. Hit
transcripts missing from the table follow ``absent``: with count the gene
counts, with ignore they add nothing to its value.
"""
import gzip
import math
import re
from pathlib import Path

import pandas as pd


FIELDS = ('table', 'level', 'gene_column', 'transcript_column', 'value_column', 'min_value',
          'min_fraction_of_target', 'absent')
REQUIRED = ('table', 'min_value')
# What the table's rows are: genes, or transcripts that the genome annotation assigns to genes.
LEVELS = ('gene', 'transcript')
# What off-target genes missing from the table do: count as expressed, or not count.
ABSENT = ('count', 'ignore')
_ENSEMBL_VERSION = re.compile(r'^(ENS[A-Z]*G\d+)\.\d+')
# A trailing version, kept apart from a GENCODE _PAR_Y suffix: ENST00000381192.10_PAR_Y.
_VERSION = re.compile(r'^(.+?)\.\d+((?:_PAR_[XY])?)$')
# Identifiers that tell what a table lists: Ensembl or RefSeq transcripts, Ensembl genes.
_TRANSCRIPT_ID = re.compile(r'^(ENS[A-Z]*T\d{6,}|[NX][MR]_\d+)')
_GENE_ID = re.compile(r'^ENS[A-Z]*G\d{6,}')


def normalize_gene_id(value) -> str:
    """Strip whitespace and an Ensembl gene version: ENSG00000075624.17 -> ENSG00000075624."""
    return _ENSEMBL_VERSION.sub(r'\1', str(value).strip())


def transcript_key(value) -> str:
    """A transcript id as listed: whitespace stripped, and only the first field of a GENCODE FASTA header."""
    return str(value).strip().split('|', 1)[0]


def unversioned(transcript_id: str) -> str:
    """ENSMUST00000073605.15 -> ENSMUST00000073605, NM_008084.3 -> NM_008084."""
    return _VERSION.sub(r'\1\2', transcript_id)


def identifier_level(values):
    """'transcript' or 'gene' when nearly all identifiers look like Ensembl/RefSeq transcript or Ensembl gene ids."""
    values = [str(v).strip() for v in values if isinstance(v, str) and v.strip()][:2000]
    if not values:
        return None
    for level, pattern in (('transcript', _TRANSCRIPT_ID), ('gene', _GENE_ID)):
        if sum(bool(pattern.match(v)) for v in values) >= 0.9 * len(values):
            return level
    return None


def validate_gene_expression(block, where='gene_expression'):
    """Check one gene_expression mapping; return it normalized (gene_columns is a list or None)."""
    if not isinstance(block, dict):
        raise ValueError(f'{where} must be a mapping with fields: {", ".join(FIELDS)}')
    unknown = [str(key) for key in block if key not in FIELDS]
    if unknown:
        raise ValueError(f'{where}: unknown field(s) {", ".join(unknown)}; expected: {", ".join(FIELDS)}')
    missing = [key for key in REQUIRED if block.get(key) in (None, '')]
    if missing:
        raise ValueError(f'{where}: missing {" and ".join(missing)}')
    table = block['table']
    if not isinstance(table, str):
        raise ValueError(f'{where}.table must be a file path, got {table!r}')
    level = block.get('level', LEVELS[0])
    if level not in LEVELS:
        raise ValueError(f'{where}.level must be one of {", ".join(LEVELS)}, got {level!r}')
    min_value = block['min_value']
    if isinstance(min_value, bool) or not isinstance(min_value, (int, float)) or not math.isfinite(min_value):
        raise ValueError(f'{where}.min_value must be a number, got {min_value!r}')
    fraction = block.get('min_fraction_of_target')
    if fraction is not None and (isinstance(fraction, bool) or not isinstance(fraction, (int, float))
                                 or not 0 <= fraction <= 1):
        raise ValueError(f'{where}.min_fraction_of_target must be a number from 0 to 1, got {fraction!r}')
    absent = block.get('absent', ABSENT[0])
    if absent not in ABSENT:
        raise ValueError(f'{where}.absent must be one of {", ".join(ABSENT)}, got {absent!r}')
    columns = block.get('gene_column')
    if columns is not None:
        if level == 'transcript':
            raise ValueError(f'{where}.gene_column applies to level: gene; with level: transcript the genome '
                             'annotation assigns transcripts to genes, and transcript_column names the id column')
        columns = [columns] if isinstance(columns, str) else columns
        if not isinstance(columns, list) or not columns or not all(isinstance(c, str) and c for c in columns):
            raise ValueError(f'{where}.gene_column must be a column name or a list of column names, '
                             f'got {block["gene_column"]!r}')
    transcript_column = block.get('transcript_column')
    if transcript_column is not None:
        if level != 'transcript':
            raise ValueError(f'{where}.transcript_column applies to level: transcript; set level: transcript '
                             'for a table of transcripts')
        if not isinstance(transcript_column, str) or not transcript_column:
            raise ValueError(f'{where}.transcript_column must be a column name, got {transcript_column!r}')
    value_column = block.get('value_column')
    if value_column is not None and (not isinstance(value_column, str) or not value_column):
        raise ValueError(f'{where}.value_column must be a column name, got {value_column!r}')
    return {'table': table, 'level': level, 'gene_columns': columns, 'transcript_column': transcript_column,
            'value_column': value_column, 'min_value': float(min_value),
            'min_fraction_of_target': None if fraction is None else float(fraction), 'absent': absent}


def _preamble_lines(path: Path) -> int:
    """Lines before the header: a GCT file's version and size lines, or leading '#' comments."""
    opener = gzip.open if path.name.lower().endswith('.gz') else open
    with opener(path, 'rt') as f:
        lines = [line for _, line in zip(range(1000), f)]
    if lines and lines[0].startswith('#1.'):
        return 2
    return next((i for i, line in enumerate(lines) if not line.startswith('#')), 0)


def read_table(path) -> pd.DataFrame:
    """Read a TSV/CSV expression table (gzipped or not) as text, skipping a GCT or '#' preamble."""
    path = Path(path)
    sep = ',' if '.csv' in path.name.lower() else '\t'
    return pd.read_csv(path, sep=sep, dtype=str, skiprows=_preamble_lines(path),
                       keep_default_na=False, na_values=[''])


def _numeric(series: pd.Series) -> bool:
    values = series.dropna()
    return len(values) > 0 and pd.to_numeric(values, errors='coerce').notna().all()


def resolve_columns(df: pd.DataFrame, settings: dict):
    """(identifier columns, value column) of a table; ValueError names what is wrong.

    The identifier columns are the gene columns, or the transcript column with level: transcript.
    """
    columns = [str(c) for c in df.columns]
    if not columns:
        raise ValueError('the expression table has no columns')
    if settings.get('level') == 'transcript':
        option = 'transcript_column'
        id_columns = [settings['transcript_column']] if settings['transcript_column'] else columns[:1]
    else:
        option = 'gene_column'
        id_columns = settings['gene_columns'] or columns[:1]
    missing = [c for c in id_columns if c not in columns]
    if missing:
        raise ValueError(f'{option} {", ".join(missing)} not in the expression table; columns: {", ".join(columns)}')
    value_column = settings['value_column']
    if value_column is None:
        numeric = [c for c in columns if c not in id_columns and _numeric(df[c])]
        if len(numeric) != 1:
            found = ', '.join(numeric) if numeric else 'none'
            raise ValueError(f'set value_column: the expression table has {len(numeric)} numeric columns ({found})')
        value_column = numeric[0]
    elif value_column not in columns:
        raise ValueError(f'value_column {value_column} not in the expression table; columns: {", ".join(columns)}')
    elif not _numeric(df[value_column]):
        raise ValueError(f'value_column {value_column} of the expression table is not numeric')
    return id_columns, value_column


class ExpressionTable:
    """Expression values of genes (by Ensembl gene id, then gene name) or, with level: transcript, of transcripts."""

    def __init__(self, block):
        self.settings = validate_gene_expression(block)
        path = Path(self.settings['table'])
        if not path.is_file():
            raise FileNotFoundError(f'gene_expression.table not found: {path}')
        df = read_table(path)
        self.level = self.settings['level']
        id_columns, self.value_column = resolve_columns(df, self.settings)
        self.gene_columns = id_columns if self.level == 'gene' else None
        self.transcript_column = id_columns[0] if self.level == 'transcript' else None
        self.min_value = self.settings['min_value']
        self.min_fraction_of_target = self.settings['min_fraction_of_target']
        self.absent = self.settings['absent']
        values = pd.to_numeric(df[self.value_column], errors='coerce')
        self.values = {}
        self.identifiers = identifier_level(df[id_columns[0]].tolist())
        normalize = normalize_gene_id if self.level == 'gene' else transcript_key
        for column in id_columns:
            for key, value in zip(df[column], values):
                if pd.isna(key) or pd.isna(value) or not str(key).strip():
                    continue
                key = normalize(key)
                # An identifier listed twice (e.g. one symbol for several ids) counts as its highest value.
                self.values[key] = max(value, self.values.get(key, value))
        # Transcripts also by id without version; None where two listed ids share it.
        self.unversioned = {}
        if self.level == 'transcript':
            for key, value in self.values.items():
                base = unversioned(key)
                self.unversioned[base] = value if base not in self.unversioned else None
        self.rows = len(df)

    def value(self, gene_id, gene_name=None):
        """Expression of a gene, or None when the table does not list it."""
        value = self.values.get(normalize_gene_id(gene_id))
        if value is None and gene_name:
            value = self.values.get(str(gene_name).strip())
        return value

    def transcript_value(self, transcript_id):
        """Expression of a transcript (level: transcript), matched with or without version; None if not listed."""
        transcript_id = transcript_key(transcript_id)
        value = self.values.get(transcript_id)
        if value is None:
            value = self.unversioned.get(unversioned(transcript_id))
        return value

    def transcripts_value(self, transcript_ids):
        """(sum over the listed transcripts or None if none is listed, [transcripts not listed])."""
        total, listed, missing = 0.0, False, []
        for transcript_id in transcript_ids:
            value = self.transcript_value(transcript_id)
            if value is None:
                missing.append(transcript_id)
            else:
                total += value
                listed = True
        return (total if listed else None), missing

    def describe(self) -> dict:
        columns = ({'gene_columns': self.gene_columns} if self.level == 'gene'
                   else {'transcript_column': self.transcript_column})
        return {'table': self.settings['table'], 'level': self.level, **columns,
                'value_column': self.value_column, 'min_value': self.min_value,
                'min_fraction_of_target': self.min_fraction_of_target, 'absent': self.absent,
                'table_rows': self.rows, 'table_identifiers': len(self.values)}

    def threshold(self, target_value) -> float:
        """Value an off-target gene of a target must reach to count; target_value None means unknown."""
        if self.min_fraction_of_target is None or target_value is None:
            return self.min_value
        return max(self.min_value, self.min_fraction_of_target * target_value)
