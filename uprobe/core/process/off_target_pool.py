"""Pool-level off-target limit: post_process.off_target_pool.

In smFISH a false spot needs many probes of a set on the same off-target RNA,
so a few probes hitting one off-target gene are harmless::

    post_process:
      off_target_pool:
        attribute: target_mapped_genes   # a mapped_genes attribute, or a list of them
        max_probes_per_gene: 2           # probes of one target's set that may hit the same gene

For each target, probes are visited from best to worst rank (the order of
post_process.sorts; the current order without sorts), and a probe is dropped
when one of its off-target genes (<attribute>_off_targets) is already hit by
max_probes_per_gene kept probes. With a gene_expression table, only expressed
off-target genes are listed and limited. ``max_probes_per_gene: 0`` keeps only
probes without off-target genes.

post_process runs it after remove_overlap and before equal_space, so the
limit applies to probes that can end up together, and equal_space chooses the
final N among the probes that respect it (a subset still does).
"""
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

from ..offtarget import OFF_TARGETS_SUFFIX, parse_gene_hits
from ..utils import get_logger

log = get_logger(__name__)

OPTIONS = ('attribute', 'max_probes_per_gene')


def pool_settings(config) -> tuple:
    """(attribute names, max_probes_per_gene) of an off_target_pool mapping; ValueError if invalid."""
    if not isinstance(config, dict):
        raise ValueError('off_target_pool must be a mapping with attribute and max_probes_per_gene')
    unknown = [str(key) for key in config if key not in OPTIONS]
    if unknown:
        raise ValueError(f'off_target_pool: unknown option(s) {", ".join(unknown)}; expected: {", ".join(OPTIONS)}')
    names = config.get('attribute')
    names = [names] if isinstance(names, str) else names
    if not isinstance(names, list) or not names or not all(isinstance(n, str) and n for n in names):
        raise ValueError(f'off_target_pool.attribute must name a mapped_genes attribute (or list several), '
                         f'got {config.get("attribute")!r}')
    limit = config.get('max_probes_per_gene')
    if type(limit) is not int or limit < 0:
        raise ValueError(f'off_target_pool.max_probes_per_gene must be an integer >= 0, got {limit!r}')
    return names, limit


def check_pool_attributes(names, attributes) -> None:
    """ValueError unless every name is a mapped_genes attribute of the protocol."""
    for name in names:
        attr = (attributes or {}).get(name)
        if not isinstance(attr, dict) or str(attr.get('type', '')).lower() != 'mapped_genes':
            raise ValueError(f'off_target_pool.attribute {name!r} is not a mapped_genes attribute')


def rank_order(df: pd.DataFrame, sorts) -> np.ndarray:
    """Row positions from best to worst by post_process.sorts (stable; current order without sorts)."""
    sorts = sorts if isinstance(sorts, dict) else {}
    keys, ascending = [], []
    for option, direction in (('is_ascending', True), ('is_descending', False)):
        for key in sorts.get(option) or []:
            if key in df.columns and key not in keys:
                keys.append(key)
                ascending.append(direction)
    positions = pd.RangeIndex(len(df))
    if not keys:
        return positions.to_numpy()
    ranked = df[keys].set_axis(positions).sort_values(keys, ascending=ascending, kind='stable',
                                                      na_position='last')
    return ranked.index.to_numpy()


def limit_off_target_pool(df: pd.DataFrame, config: dict, sorts=None):
    """Drop the lower-ranked probes of each target beyond max_probes_per_gene per off-target gene.

    Returns (kept rows in their original order, report for summary.json).
    """
    names, limit = pool_settings(config)
    columns = [name + OFF_TARGETS_SUFFIX for name in names]
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f'off_target_pool needs {", ".join(missing)}, written by mapped_genes attributes')
    target_column = next((c for c in ('target', 'gene') if c in df.columns), None)
    targets = df[target_column].astype(str).to_numpy() if target_column else np.full(len(df), '')
    genes = [set().union(*(parse_gene_hits(cell) for cell in cells)) for cells in zip(*(df[c] for c in columns))]
    hits = defaultdict(int)  # (target, gene) -> kept probes that hit it
    keep = np.zeros(len(df), dtype=bool)
    blocked = Counter()
    for position in rank_order(df, sorts):
        target = targets[position]
        full = [gene for gene in genes[position] if hits[(target, gene)] >= limit]
        if full:
            blocked.update(full)
            continue
        keep[position] = True
        for gene in genes[position]:
            hits[(target, gene)] += 1
    removed = int((~keep).sum())
    report = {'attributes': names, 'max_probes_per_gene': limit, 'removed': removed,
              'genes_at_limit': dict(blocked.most_common(20))}
    if removed:
        top = ', '.join(f'{gene} ({n})' for gene, n in blocked.most_common(5))
        log.info(f'off_target_pool: removed {removed} probes whose off-target genes were already hit by '
                 f'{limit} better-ranked probes of the same target; most often {top}')
    return df.iloc[keep], report
