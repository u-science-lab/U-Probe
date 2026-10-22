import pandas as pd
import numpy as np
from typing import Dict, Any, Optional
from uprobe.core.utils import get_logger
from .regions import region_bounds
from . import oligo_sites

log = get_logger(__name__)


def number_desired(config: Dict[str, Any], target: str) -> Optional[int]:
    """Probes to keep for one target: its own entry, else the global number_desired (None: keep all)."""
    if not isinstance(config, dict):
        return None
    entry = config.get(target)
    if isinstance(entry, dict) and 'number_desired' in entry:
        return entry['number_desired']
    return config.get('number_desired')


def positions_along_target(df: pd.DataFrame) -> pd.Series:
    """Start of each candidate along its target.

    mRNA, exon, CDS and UTR windows number start/end within their exon or
    UTR (a window across exon junctions: along the exons it spans), so they
    carry target_start, their position on the gene's merged exon model, or on
    its one selected transcript (gen.fun.SPLICED_SOURCES).
    Genome windows and direct sequences are numbered along the region or
    sequence itself.
    """
    if 'target_start' in df.columns:
        start = pd.to_numeric(df['target_start'], errors='coerce')
        return start.fillna(region_bounds(df)['start']) if start.isna().any() else start
    names = [df[c] for c in ('exon_name', 'utr_name') if c in df.columns]
    if names and pd.concat(names).dropna().astype(str).ne('.').any():
        log.warning("equal_space: no target_start column (a table from an older U-Probe?); exon and UTR "
                    "windows are ordered by their position within the exon or UTR")
    return region_bounds(df)['start']


def _collapse_identical_oligos(target_df: pd.DataFrame):
    """One row per distinct oligo of the target (the first site along it), so that an oligo that
    binds several sites counts once; a no-op without duplicate window sequences."""
    work = target_df.reset_index(drop=True)
    ids = oligo_sites.oligo_ids(work, ['target'])
    if ids is None or not ids.duplicated().any():
        return target_df, []
    bounds = pd.DataFrame({'start': positions_along_target(work)})
    bounds['end'] = bounds['start'] + work[oligo_sites.OLIGO_COLUMN].str.len() - 1
    result, report = oligo_sites.collapse_identical(work, ['target'], bounds)
    log.info(f"equal_space: {len(work) - len(result)} candidate(s) repeat another one's sequence on the same "
             f"target and were counted once")
    return result.reset_index(drop=True), report


def equal_space(df: pd.DataFrame, config: Dict[str, Any]) -> pd.DataFrame:
    """Keep number_desired probes per target, spread evenly along the target.

    Candidates are ordered along the target and every k-th one is kept
    (evenly spaced ranks, both ends included). A target with fewer passing
    candidates keeps all of them.
    """
    if df.empty:
        return df
    if 'target' not in df.columns:
        log.warning("DataFrame does not have 'target' column, skip equal_space processing")
        return df
    result_dfs = []
    multi_site = []
    targets = df['target'].unique()
    log.info(f"equal_space processing started, there are {len(targets)} targets")
    for target in targets:
        target_df = df[df['target'] == target].copy()
        target_df, multi = _collapse_identical_oligos(target_df)
        multi_site += multi
        current_count = len(target_df)
        desired = number_desired(config, target)
        if desired is None:
            log.warning(f"target {target} does not have number_desired, skip processing")
            result_dfs.append(target_df)
            continue
        
        if current_count <= desired:
            if current_count < desired:
                log.warning(f"target {target}: only {current_count} probes passed, fewer than the "
                            f"{desired} requested; keeping all of them")
            else:
                log.info(f"target {target}: keeping all {current_count} probes (desired: {desired})")
            result_dfs.append(target_df)
        else:
            # sort by location along the target and select evenly
            if any(column in target_df.columns for column in ('target_start', 'sub_region', 'start')):
                target_df['_start_pos'] = positions_along_target(target_df)
                target_df = target_df.sort_values('_start_pos', kind='stable')
                target_df = target_df.drop('_start_pos', axis=1)
            indices = np.linspace(0, current_count - 1, desired, dtype=int)
            selected_df = target_df.iloc[indices].copy()
            log.info(f"target {target}: selected {len(selected_df)} probes from {current_count} (desired: {desired})")
            result_dfs.append(selected_df)
    if result_dfs:
        final_df = pd.concat(result_dfs, ignore_index=True)
        if multi_site:
            final_df.attrs['multi_site_oligos'] = multi_site
        log.info(f"equal_space processing completed, final number of probes: {len(final_df)}")
        return final_df
    else:
        return df
