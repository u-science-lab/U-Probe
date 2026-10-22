"""Identical oligos that bind several sites of one target (internal repeats).

A target with an internal repeat (tandem copies of a domain, a repeated
exon) yields several candidate windows with the same sequence. They are one
oligo: it is counted once, and every site it binds is occupied, so
remove_overlap and equal_space treat the sites of an oligo together.
Oligos with a single site are unaffected.

The kept row of an oligo with several sites carries ``n_sites`` and
``oligo_sites``; ``attrs['multi_site_oligos']`` and the run summary list the
oligos of the finished table that do (see settle), each once.
"""
from bisect import bisect_right

import pandas as pd

OLIGO_COLUMN = 'target_region'  # the target window, from which the binding sequence derives
N_SITES = 'n_sites'
SITES = 'oligo_sites'
REPORT = 'multi_site_oligos'


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
    keys = pd.DataFrame({i: column.to_numpy() for i, column in enumerate(columns)}, index=df.index)
    return keys.groupby(list(keys.columns), sort=False, dropna=False).ngroup()


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
    """Accepted sites of one target, for conflict checks in O(log n) each.

    Accepted sites never overlap (a site is accepted only when it is more than
    the gap away from every accepted one, and the sites of one oligo are
    merged first), so sorted by start they are sorted by end too, and the last
    site starting at or before a candidate's end + gap is the only one to check.
    """

    def __init__(self, location_interval: int):
        self.gap = location_interval
        self.starts, self.ends = [], []

    def conflicts(self, sites) -> bool:
        starts, ends, gap = self.starts, self.ends, self.gap
        for start, end in sites:
            i = bisect_right(starts, end + gap) - 1  # last site starting early enough
            if i >= 0 and ends[i] + gap >= start:
                return True
        return False

    def add(self, sites):
        for start, end in sites:
            i = bisect_right(self.starts, start)
            self.starts.insert(i, start)
            self.ends.insert(i, end)


def sites_by_oligo(ids, bounds: pd.DataFrame) -> dict:
    """{oligo id: its distinct (start, end) sites, sorted}."""
    sites_of = {}
    for oligo, start, end in zip(list(ids), bounds['start'].tolist(), bounds['end'].tolist()):
        site = (int(start), int(end))
        found = sites_of.get(oligo)
        if found is None:
            sites_of[oligo] = site  # one site, most oligos
        elif type(found) is tuple:
            if found != site:
                sites_of[oligo] = {found, site}
        else:
            found.add(site)
    return {oligo: [sites] if type(sites) is tuple else sorted(sites) for oligo, sites in sites_of.items()}


def annotate(result: pd.DataFrame, kept, ids, sites_of) -> pd.DataFrame:
    """Add n_sites / oligo_sites to the kept rows (result is aligned with kept)."""
    result[N_SITES] = [len(sites_of[ids[label]]) for label in kept]
    result[SITES] = [';'.join(f'{s}-{e}' for s, e in sites_of[ids[label]]) for label in kept]
    return result


def annotate_kept(result: pd.DataFrame, kept, ids, sites_of) -> pd.DataFrame:
    """annotate, only when a kept oligo binds several sites (otherwise the table is unchanged)."""
    if any(len(sites_of[ids[label]]) > 1 for label in kept):
        return annotate(result, kept, ids, sites_of)
    return result


def report_from_table(df: pd.DataFrame) -> list:
    """[{target, sequence, n_sites, sites}] of the rows of ``df`` whose oligo binds several sites, each oligo once."""
    if N_SITES not in df.columns or SITES not in df.columns or OLIGO_COLUMN not in df.columns:
        return []
    report, index = [], {}
    counts = pd.to_numeric(df[N_SITES], errors='coerce')
    targets = df['target'] if 'target' in df.columns else pd.Series(None, index=df.index)
    for count, sites, sequence, target in zip(counts, df[SITES], df[OLIGO_COLUMN], targets):
        if not count > 1 or not isinstance(sites, str):
            continue
        target = None if target is None or pd.isna(target) else str(target)
        key = (target, sequence)
        if key in index:  # the same oligo twice (layouts, or groups of an older table): one entry, all sites
            entry = index[key]
            entry['sites'] = [f'{s}-{e}' for s, e in sorted(set(map(_site, entry['sites'] + sites.split(';'))))]
            entry['n_sites'] = len(entry['sites'])
            continue
        index[key] = {'target': target, 'sequence': sequence, 'n_sites': int(count), 'sites': sites.split(';')}
        report.append(index[key])
    return report


def _site(text: str) -> tuple:
    start, end = text.split('-')
    return int(start), int(end)


def settle(df: pd.DataFrame) -> pd.DataFrame:
    """Report the multi-site oligos of a finished table in attrs['multi_site_oligos'] (from its own rows);
    drop n_sites and oligo_sites when none of its rows binds several sites."""
    df.attrs.pop(REPORT, None)
    report = report_from_table(df)
    if report:
        df.attrs[REPORT] = report
    elif N_SITES in df.columns or SITES in df.columns:
        df = df.drop(columns=[c for c in (N_SITES, SITES) if c in df.columns])
        df.attrs.pop(REPORT, None)
    return df


def collapse_identical(df: pd.DataFrame, group_cols: list, bounds: pd.DataFrame, transcript_col: str = None):
    """One row per oligo (the first site along the target), annotated when an oligo has several sites.

    With ``transcript_col``, the kept row lists the transcripts of every row of its oligo at its site
    (windows that several transcripts or exon entries share at one position), as remove_overlap does.
    """
    ids = oligo_ids(df, group_cols)
    if ids is None or not ids.duplicated().any():
        return df
    order = bounds['start'].sort_values(kind='stable').index
    first = set(ids.loc[order].drop_duplicates().index)
    kept = [label for label in df.index if label in first]
    sites_of = sites_by_oligo(ids, bounds)
    result = annotate_kept(df.loc[kept].copy(), kept, ids, sites_of)
    if transcript_col:
        merge_transcripts(df, result, kept, ids, bounds, transcript_col)
    return result


def _names(value) -> list:
    if isinstance(value, (list, tuple)):
        return [str(name) for name in value]
    if isinstance(value, str):
        return [name.strip() for name in value.split(',') if name.strip()]
    return []


def merge_transcripts(df, result, kept, ids, bounds, transcript_col):
    """List on each kept row the transcripts of all rows of its oligo at its site (in place)."""
    names = {}
    for label, oligo, start, end, cell in zip(df.index, ids, bounds['start'], bounds['end'], df[transcript_col]):
        names.setdefault((oligo, start, end), []).extend(_names(cell))
    for label in kept:
        cell = result.at[label, transcript_col]
        own = _names(cell)
        merged = list(dict.fromkeys(own + names[(ids[label], bounds.at[label, 'start'], bounds.at[label, 'end'])]))
        if len(merged) > len(own):
            result.at[label, transcript_col] = merged if isinstance(cell, (list, tuple)) else ', '.join(merged)
            if 'n_trans' in result.columns:
                result.at[label, 'n_trans'] = len(merged)
