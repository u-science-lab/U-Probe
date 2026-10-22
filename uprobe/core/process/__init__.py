import pandas as pd 
import re
from .otp import avoid_otp
from .equal_space import equal_space
from .off_target_pool import limit_off_target_pool
from .regions import region_bounds
from .summary import process_summary
from . import oligo_sites
from .sequence_filter import apply_sequence_filter
from uprobe.core.utils import get_logger

logger = get_logger(__name__)


def parse_condition(condition: str) -> list:
    """Parse a condition string into list of (column, operator, value) tuples
    """
    conditions = [c.strip() for c in condition.split('&')]
    result = []
    # Match patterns like "column <= 80" or "column >= 35"
    pattern = r'(\w+)\s*([<>=]+)\s*([+-]?\d+(?:\.\d+)?)'
    for cond in conditions:
        match = re.match(pattern, cond)
        if not match:
            raise ValueError(f"Invalid condition format: {cond}")
        column, operator, value = match.groups()
        try:
            value = int(value)
        except ValueError:
            value = float(value)
        result.append((column, operator, value))
    return result

def apply_condition(df: pd.DataFrame,
                    conditions: list
                    ) -> pd.DataFrame:
    for column, operator, value in conditions:
        if column not in df.columns:
            continue
        if operator == '<=':
            df = df[df[column] <= value]
        elif operator == '>=':
            df = df[df[column] >= value]
        elif operator == '<':
            df = df[df[column] < value]
        elif operator == '>':
            df = df[df[column] > value]
        elif operator == '==':
            df = df[df[column] == value]
        else:
            raise ValueError(f"Unsupported operator: {operator}")
    return df

def filter_table(df: pd.DataFrame,
                  filters: dict
                  ) -> pd.DataFrame:
    for filter_name, filter_config in filters.items():
        if isinstance(filter_config, dict) and filter_config.get('type') == 'sequence_pattern':
            df = apply_sequence_filter(df, filter_config)
            continue
        if isinstance(filter_config, dict) and 'condition' in filter_config:
            condition = filter_config['condition']
            if isinstance(condition, str):
                conditions = parse_condition(condition)
                df = apply_condition(df, conditions)
            elif isinstance(condition, list):
                for cond in condition:
                    conditions = parse_condition(cond)
                    df = apply_condition(df, conditions)
    return df

def sort_table(df: pd.DataFrame, 
               keys: list, 
               ascending: list
               ) -> pd.DataFrame:
    valid_keys = [key for key in keys if key in df.columns]
    valid_ascending = ascending[:len(valid_keys)]
    if not valid_keys:
        return df  
    return df.sort_values(by=valid_keys, ascending=valid_ascending)

def _split_names(value) -> list:
    """Transcript names of a 'T1, T2' cell."""
    if not isinstance(value, str):
        return []
    return [name.strip() for name in value.split(',') if name.strip()]


def _non_overlapping(df: pd.DataFrame, group_cols: list, bounds: pd.DataFrame,
                     location_interval: int, prefer: str = None, ids=None) -> list:
    """Index labels of the rows to keep, group by group.

    Rows are visited by start, and a row is kept when it starts more than
    location_interval bases after the end of the last kept row. Ties at a
    start (the layouts of a window with length ranges) are visited in table
    order, where the shorter windows come first, or, with ``prefer``, by
    that column ascending, so the start keeps its best candidate.

    ``ids`` (oligo_sites.oligo_ids) names the identical oligos: one that binds
    several sites is kept once, only when none of its sites is within
    location_interval of a kept site, and then occupies all of them.
    """
    work = df[group_cols].assign(_start=bounds['start'], _end=bounds['end'])
    order = group_cols + ['_start']
    if prefer:
        work['_prefer'] = pd.to_numeric(df[prefer], errors='coerce')
        order.append('_prefer')
    work = work.sort_values(by=order, kind='stable', na_position='last')
    if ids is not None and not ids.duplicated().any():
        ids = None  # every oligo has one site
    keep = []
    for _, group in work.groupby(group_cols, dropna=False, sort=False):
        if ids is None:
            current_end = None
            for label, start, end in zip(group.index, group['_start'], group['_end']):
                if current_end is None or start > current_end + location_interval:
                    keep.append(label)
                    current_end = end
            continue
        sites_of = oligo_sites.sites_by_oligo(ids.loc[group.index], bounds.loc[group.index])
        occupied = oligo_sites.Occupancy(location_interval)
        decided = set()
        for label in group.index:
            oligo = ids[label]
            if oligo in decided:
                continue
            decided.add(oligo)  # refused once, refused later: the occupied sites only grow
            sites = oligo_sites.merged_sites(sites_of[oligo])
            if not occupied.conflicts(sites):
                keep.append(label)
                occupied.add(sites)
    return keep


def _with_site_report(df, result, keep, ids, bounds):
    """Name the kept oligos that bind several sites: n_sites and oligo_sites columns (only then)
    and result.attrs['multi_site_oligos']."""
    if ids is None or not ids.duplicated().any():
        return result
    report, sites_of = oligo_sites.describe_kept(df, ids, bounds, keep)
    if report:
        result = oligo_sites.annotate(result, keep, ids, sites_of)
        result.attrs['multi_site_oligos'] = report
    return result


def _remove_overlap_along_targets(df: pd.DataFrame, transcript_col, location_interval: int,
                                  prefer: str = None) -> pd.DataFrame:
    """remove_overlap for mRNA, exon, CDS and UTR windows, on each target's merged exon model.

    target_start/target_end place every window on the gene's merged exon
    model (gen.fun.segment_offsets), or on its one selected transcript,
    junction windows included (gen.fun.SPLICED_SOURCES), so windows of
    alternative exon entries of different transcripts are compared, and a
    site that several entries contain is kept once; the kept row lists the
    transcripts of all of them.
    """
    bounds = pd.DataFrame({'start': pd.to_numeric(df['target_start'], errors='coerce'),
                           'end': pd.to_numeric(df['target_end'], errors='coerce')}, index=df.index)
    if bounds.isna().any().any():  # direct sequences: positions along the sequence itself
        bounds = bounds.fillna(region_bounds(df))
    ids = oligo_sites.oligo_ids(df, ['target'])
    keep = _non_overlapping(df, ['target'], bounds, location_interval, prefer, ids)
    result = _with_site_report(df, df.loc[keep].copy(), keep, ids, bounds)
    if transcript_col:
        sites = list(zip(df['target'], bounds['start'], bounds['end']))
        names = {}
        for site, cell in zip(sites, df[transcript_col]):
            names.setdefault(site, set()).update(_split_names(cell))
        site_of = dict(zip(df.index, sites))
        for label in keep:
            merged = sorted(names[site_of[label]])
            if len(merged) > len(_split_names(result.at[label, transcript_col])):
                result.at[label, transcript_col] = ', '.join(merged)
                if 'n_trans' in result.columns:
                    result.at[label, 'n_trans'] = len(merged)
    return result.reset_index(drop=True)


def remove_overlap(df: pd.DataFrame, 
                   config: dict,
                   location_interval: int,
                   prefer: str = None,
                   ) -> pd.DataFrame:
    """Keep windows that do not overlap, along each target; ``prefer`` (a column, lower is
    better) picks among the candidates of one start, which keep the shortest otherwise."""
    if prefer and prefer not in df.columns:
        raise ValueError(f"post_process.remove_overlap.prefer: no column {prefer!r}")
    df = df.reset_index(drop=True)
    transcript_col = None
    if 'transcript_name' in df.columns:
        transcript_col = 'transcript_name'
    elif 'transcript_names' in df.columns:
        transcript_col = 'transcript_names'
    if transcript_col:
        df[transcript_col] = df[transcript_col].apply(lambda x: 
                                                     ', '.join(x) if isinstance(x, list) else x)
    if 'target_start' in df.columns and 'target' in df.columns:
        return _remove_overlap_along_targets(df, transcript_col, location_interval, prefer)
    if transcript_col:
        # Positions are per target and per extracted exon or UTR: exon/UTR
        # windows are numbered within their exon or UTR, and direct-sequence
        # targets all share the transcript name 'custom'. Keep rows whose
        # exon_name or utr_name is missing (genome and direct targets mixed).
        region_cols = [col for col in ('exon_name', 'utr_name') if col in df.columns]
        group_cols = (['target'] if 'target' in df.columns else []) + region_cols + [transcript_col]
    elif 'target' in df.columns:
        group_cols = ['target']
    else:
        return df
    # Numeric bounds: genome windows carry only sub_region ("101-140"), also
    # next to direct sequences, whose start/end columns they leave empty.
    bounds = region_bounds(df)
    ids = oligo_sites.oligo_ids(df, group_cols)
    keep = _non_overlapping(df, group_cols, bounds, location_interval, prefer, ids)
    return _with_site_report(df, df.loc[keep], keep, ids, bounds).reset_index(drop=True)

# post_process runs its steps in this order, whatever their order in the protocol.
STEP_ORDER = ('filters', 'avoid_otp', 'remove_overlap', 'off_target_pool', 'equal_space', 'sorts')

def post_process(df: pd.DataFrame, 
                 config: dict
                 ) -> pd.DataFrame:
    """Apply the post_process steps in STEP_ORDER.

    off_target_pool limits the probes per off-target gene among the
    non-overlapping candidates, before equal_space keeps N per target.
    """
    stages = []
    pool_report = None
    multi_site = []  # identical oligos that bind several sites of a target (oligo_sites.py)
    # Per-target counts after each stage (aligned with `stages`) for the run summary.
    target_counts = []
    target_col = next((c for c in ('target', 'gene') if c in df.columns), None)
    previous_count = len(df)
    def record(name):
        nonlocal previous_count
        stages.append({'stage': name, 'before': previous_count, 'after': len(df), 'removed': previous_count - len(df)})
        if target_col:
            counts = df[target_col].dropna().astype(str).value_counts() if target_col in df.columns else {}
            target_counts.append({str(k): int(v) for k, v in counts.items()})
        previous_count = len(df)
    processes = config.get('post_process', {})
    if 'filters' in processes and processes['filters']:
        logger.info("Filtering the table")
        filters = processes['filters']
        df = filter_table(df, filters)
        record('Filters')
    if "avoid_otp" in processes and processes['avoid_otp']:
        logger.info("Avoiding OTP")
        config = processes['avoid_otp']
        mapped_sites_cols = []
        for col in df.columns:
            if not col.endswith('_num') and f"{col}_num" in df.columns:
                sample_val = None
                non_na_series = df[col].dropna()
                for val in non_na_series:
                    try:
                        if isinstance(val, str):
                            import ast
                            parsed_val = ast.literal_eval(val)
                            if isinstance(parsed_val, list) and len(parsed_val) > 0:
                                sample_val = parsed_val
                                break
                        elif isinstance(val, list) and len(val) > 0:
                            sample_val = val
                            break
                    except Exception:
                        continue
                if sample_val is not None:
                    if isinstance(sample_val[0], (tuple, list)) and len(sample_val[0]) == 3:
                        mapped_sites_cols.append(col)
        if mapped_sites_cols:
            for target, target_config in config.items():
                target_regions = target_config.get('target_regions', [])
                density_thresh = float(target_config.get('density_thresh', 1e-4))
                search_range = tuple(map(float, target_config.get('search_range', [-1000000, 1000000])))
                avoid_target_overlap = target_config.get('avoid_target_overlap', True)
                target_df = df[df['target'] == target] if 'target' in df.columns else df
                if target_df.empty:
                    logger.warning(f"No probes found for target {target}, skipping OTP filtering")
                    continue
                blocks = []
                for _, row in target_df.iterrows():
                    probe_id = row.get('probe_id', f"probe_{row.name}")
                    seq = row.get('target_region', '')
                    all_alns = []
                    for col in mapped_sites_cols:
                        mapped_sites = row[col]
                        try:
                            if isinstance(mapped_sites, str):
                                import ast
                                mapped_sites = ast.literal_eval(mapped_sites)
                            if isinstance(mapped_sites, list) and len(mapped_sites) > 0:
                                alns = [(ref_name, start_pos-1, end_pos) for ref_name, start_pos, end_pos in mapped_sites]
                                all_alns.extend(alns)
                        except Exception:
                            pass
                    if all_alns:
                        blocks.append((probe_id, seq, all_alns))
                if not blocks:
                    logger.warning(f"No alignment data found for target {target}, skipping OTP filtering")
                    continue
                counted = avoid_otp(blocks, target_regions, density_thresh, avoid_target_overlap, search_range)
                if counted:
                    filtered_probe_ids = []
                    on_target_counts = {}
                    off_target_counts = {}
                    on_target_rates = {}
                    for b, c in counted:
                        name, seq, alns = b
                        on_target = c[0]
                        off_target = c[1]
                        total = on_target + off_target
                        on_target_rate = on_target / total if total > 0 else 0
                        filtered_probe_ids.append(name)
                        on_target_counts[name] = on_target
                        off_target_counts[name] = off_target
                        on_target_rates[name] = on_target_rate
                    if 'probe_id' in df.columns:
                        target_mask = df['target'] == target if 'target' in df.columns else pd.Series([True] * len(df))
                        probe_mask = df['probe_id'].isin(filtered_probe_ids)
                        keep_mask = (~target_mask) | probe_mask
                        df = df[keep_mask].reset_index(drop=True)
                        for probe_id in filtered_probe_ids:
                            mask = df['probe_id'] == probe_id
                            df.loc[mask, 'on_target'] = on_target_counts[probe_id]
                            df.loc[mask, 'off_target'] = off_target_counts[probe_id]
                            df.loc[mask, 'target_ratio'] = on_target_rates[probe_id]
                else:
                    if 'target' in df.columns:
                        df = df[df['target'] != target].reset_index(drop=True)
        else:
            logger.error("No mapped_sites columns found in DataFrame, please check the input data")
        record('Avoid off-target priming')
    if "remove_overlap" in processes and processes['remove_overlap']:
        logger.info("Removing overlap")
        location_interval = processes['remove_overlap'].get('location_interval', 0)
        df = remove_overlap(df, config, location_interval, processes['remove_overlap'].get('prefer'))
        multi_site += df.attrs.get('multi_site_oligos', [])
        record('Remove overlap')
    if processes.get('off_target_pool'):
        logger.info("Limiting probes per off-target gene")
        df, pool_report = limit_off_target_pool(df, processes['off_target_pool'], processes.get('sorts'))
        record('Off-target pool limit')
    if "equal_space" in processes and processes['equal_space']:
        logger.info("Equalizing space")
        df = equal_space(df, processes['equal_space'])
        multi_site += df.attrs.get('multi_site_oligos', [])
        record('Equal spacing')
    if 'sorts' in processes and processes['sorts']:
        logger.info("Sorting the table")
        sorts = processes['sorts']
        pos_fields = sorts.get('is_ascending', [])
        neg_fields = sorts.get('is_descending', [])
        sort_keys = pos_fields + neg_fields
        is_ascending = [True] * len(pos_fields) + [False] * len(neg_fields)
        df = sort_table(df, sort_keys, is_ascending)
        record('Sort')
    if 'summary' in processes and processes['summary']:
        logger.info("Summary configuration found - will be processed during report generation")
    df.attrs['post_process_stages'] = stages
    df.attrs['post_process_target_counts'] = target_counts
    if pool_report is not None:
        df.attrs['off_target_pool'] = pool_report
    if multi_site:
        df.attrs['multi_site_oligos'] = multi_site
    return df
