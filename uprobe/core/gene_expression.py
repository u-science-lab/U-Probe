"""Expression-aware off-target counting.

A protocol may supply a gene expression table for the sample the probes are
used in::

    gene_expression:
      table: hek293_tpm.tsv    # TSV, CSV (name contains .csv) or GTEx GCT; may be gzipped
      gene_column: gene_id     # default: the first column; a list of columns is also accepted
      value_column: TPM        # default: the only numeric column
      min_value: 1             # genes with a lower value are not off-targets

``mapped_genes`` attributes then ignore off-target genes whose value is below
``min_value``: an RNA that is not there cannot take up a probe. Genes are
looked up by Ensembl gene id (with or without version) and then by gene name
(symbol), so a table keyed by either works. Off-target genes absent from the
table still count (conservative); the run reports how many there were. The
target gene always counts. Relative paths are resolved from the directory
where uprobe runs, as for genomes.yaml.
"""
import gzip
import math
import re
from pathlib import Path

import pandas as pd


FIELDS = ('table', 'gene_column', 'value_column', 'min_value')
REQUIRED = ('table', 'min_value')
_ENSEMBL_VERSION = re.compile(r'^(ENS[A-Z]*G\d+)\.\d+')


def normalize_gene_id(value) -> str:
    """Strip whitespace and an Ensembl gene version: ENSG00000075624.17 -> ENSG00000075624."""
    return _ENSEMBL_VERSION.sub(r'\1', str(value).strip())


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
    min_value = block['min_value']
    if isinstance(min_value, bool) or not isinstance(min_value, (int, float)) or not math.isfinite(min_value):
        raise ValueError(f'{where}.min_value must be a number, got {min_value!r}')
    columns = block.get('gene_column')
    if columns is not None:
        columns = [columns] if isinstance(columns, str) else columns
        if not isinstance(columns, list) or not columns or not all(isinstance(c, str) and c for c in columns):
            raise ValueError(f'{where}.gene_column must be a column name or a list of column names, '
                             f'got {block["gene_column"]!r}')
    value_column = block.get('value_column')
    if value_column is not None and (not isinstance(value_column, str) or not value_column):
        raise ValueError(f'{where}.value_column must be a column name, got {value_column!r}')
    return {'table': table, 'gene_columns': columns, 'value_column': value_column, 'min_value': float(min_value)}


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
    """(gene columns, value column) of a table; ValueError names what is wrong."""
    columns = [str(c) for c in df.columns]
    if not columns:
        raise ValueError('the expression table has no columns')
    gene_columns = settings['gene_columns'] or columns[:1]
    missing = [c for c in gene_columns if c not in columns]
    if missing:
        raise ValueError(f'gene_column {", ".join(missing)} not in the expression table; columns: {", ".join(columns)}')
    value_column = settings['value_column']
    if value_column is None:
        numeric = [c for c in columns if c not in gene_columns and _numeric(df[c])]
        if len(numeric) != 1:
            found = ', '.join(numeric) if numeric else 'none'
            raise ValueError(f'set value_column: the expression table has {len(numeric)} numeric columns ({found})')
        value_column = numeric[0]
    elif value_column not in columns:
        raise ValueError(f'value_column {value_column} not in the expression table; columns: {", ".join(columns)}')
    elif not _numeric(df[value_column]):
        raise ValueError(f'value_column {value_column} of the expression table is not numeric')
    return gene_columns, value_column


class ExpressionTable:
    """Gene -> expression value, looked up by Ensembl gene id and then by gene name."""

    def __init__(self, block):
        self.settings = validate_gene_expression(block)
        path = Path(self.settings['table'])
        if not path.is_file():
            raise FileNotFoundError(f'gene_expression.table not found: {path}')
        df = read_table(path)
        self.gene_columns, self.value_column = resolve_columns(df, self.settings)
        self.min_value = self.settings['min_value']
        values = pd.to_numeric(df[self.value_column], errors='coerce')
        self.values = {}
        for column in self.gene_columns:
            for key, value in zip(df[column], values):
                if pd.isna(key) or pd.isna(value) or not str(key).strip():
                    continue
                key = normalize_gene_id(key)
                # A gene listed twice (e.g. one symbol for several ids) counts as its highest value.
                self.values[key] = max(value, self.values.get(key, value))
        self.rows = len(df)

    def value(self, gene_id, gene_name=None):
        """Expression of a gene, or None when the table does not list it."""
        value = self.values.get(normalize_gene_id(gene_id))
        if value is None and gene_name:
            value = self.values.get(str(gene_name).strip())
        return value

    def describe(self) -> dict:
        return {'table': self.settings['table'], 'gene_columns': self.gene_columns,
                'value_column': self.value_column, 'min_value': self.min_value,
                'table_rows': self.rows, 'table_identifiers': len(self.values)}
