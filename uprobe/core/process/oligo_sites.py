"""Identical oligos that bind several sites of one target (internal repeats).

A target with an internal repeat (tandem copies of a domain, a repeated
exon) yields several candidate windows with the same sequence. They are one
oligo: it is counted once, and every site it binds is occupied, so
remove_overlap and equal_space treat the sites of an oligo together.
Oligos with a single site are unaffected.
"""
from bisect import bisect_right, insort

import pandas as pd

OLIGO_COLUMN = 'target_region'  # the target window, from which the binding sequence derives
N_SITES = 'n_sites'
SITES = 'oligo_sites'


def oligo_ids(df: pd.DataFrame, group_cols: list):
    """One id per row: rows of a group with the same window sequence (and layout) share it. None: not available."""
    if OLIGO_COLUMN not in df.columns:
        return None
    seq = df[OLIGO_COLUMN]
    if not seq.map(lambda v: isinstance(v, str) and bool(v)).all():
        return None
    # A layout (named parts of the window) is part of the oligo: two layouts of one window differ.
    columns = [df[c].astype(str) for c in group_cols] + [seq]
    if 'layout_id' in df.columns:
        columns.append(df['layout_id'].astype(str))
    seen = {}
    return pd.Series([seen.setdefault(key, len(seen)) for key in zip(*columns)], index=df.index)


def merged_sites(sites) -> list:
    """Distinct (start, end) sites, with overlapping ones of the same oligo joined."""
    merged = []
    for start, end in sorted(set(sites)):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


class Occupancy:
    """Accepted sites of one target, mutually separated, for conflict checks."""

    def __init__(self, location_interval: int):
        self.gap = location_interval
        self.sites = []  # sorted (start, end); the ends are sorted too

    def conflicts(self, sites) -> bool:
        starts = [s for s, _ in self.sites]
        for start, end in sites:
            i = bisect_right(starts, end + self.gap) - 1  # last site starting early enough
            if i >= 0 and self.sites[i][1] + self.gap >= start:
                return True
        return False

    def add(self, sites):
        for site in sites:
            insort(self.sites, site)


def sites_by_oligo(ids, bounds: pd.DataFrame) -> dict:
    sites_of = {}
    for oligo, start, end in zip(ids, bounds['start'], bounds['end']):
        sites_of.setdefault(oligo, set()).add((int(start), int(end)))
    return {oligo: sorted(sites) for oligo, sites in sites_of.items()}


def describe_kept(df: pd.DataFrame, ids, bounds: pd.DataFrame, kept) -> tuple:
    """Multi-site oligos among the kept row labels: [{target, sequence, n_sites, sites}]."""
    sites_of = sites_by_oligo(ids, bounds)
    report = []
    for label in kept:
        sites = sites_of[ids[label]]
        if len(sites) > 1:
            target = df.at[label, 'target'] if 'target' in df.columns else None
            report.append({'target': None if target is None else str(target),
                           'sequence': df.at[label, OLIGO_COLUMN], 'n_sites': len(sites),
                           'sites': [f'{s}-{e}' for s, e in sites]})
    return report, sites_of


def annotate(result: pd.DataFrame, kept, ids, sites_of) -> pd.DataFrame:
    """Add n_sites / oligo_sites to the kept rows (result is aligned with kept)."""
    result[N_SITES] = [len(sites_of[ids[label]]) for label in kept]
    result[SITES] = [';'.join(f'{s}-{e}' for s, e in sites_of[ids[label]]) for label in kept]
    return result


def collapse_identical(df: pd.DataFrame, group_cols: list, bounds: pd.DataFrame):
    """One row per oligo (the first site along the target); return (table, multi-site report)."""
    ids = oligo_ids(df, group_cols)
    if ids is None or not ids.duplicated().any():
        return df, []
    order = bounds['start'].sort_values(kind='stable').index
    first = set(ids.loc[order].drop_duplicates().index)
    kept = [label for label in df.index if label in first]
    report, sites_of = describe_kept(df, ids, bounds, kept)
    result = df.loc[kept].copy()
    if report:
        result = annotate(result, kept, ids, sites_of)
    return result, report
