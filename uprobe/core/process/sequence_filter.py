"""Configurable exclusions on any generated sequence column."""
import re


def apply_sequence_filter(df, config):
    target = config.get('target')
    if not isinstance(target, str) or not target:
        raise ValueError('Sequence filter requires a target column')
    column = target if target in df.columns else target.replace(':', '.')
    if column not in df.columns:
        raise ValueError(f'Sequence filter target column is missing: {target}')
    patterns = config.get('exclude_patterns')
    if not isinstance(patterns, list) or not patterns or any(
        not isinstance(p, str) or not p for p in patterns
    ):
        raise ValueError('Sequence filter requires non-empty exclude_patterns')
    try:
        compiled = [re.compile(p, re.IGNORECASE) for p in patterns]
    except re.error as exc:
        raise ValueError(f'Invalid sequence exclusion pattern: {exc}') from exc
    # Unavailable or empty sequences cannot pass a sequence requirement.
    valid = df[column].map(lambda s: isinstance(s, str) and bool(s.strip()))
    sequences = df[column].where(valid, '').str.strip()
    for pattern in compiled:
        valid &= ~sequences.str.contains(pattern, na=False)
    return df.loc[valid]


def sequence_filter_description(config):
    patterns = ', '.join(config.get('exclude_patterns', []))
    return f"{config.get('target', '')}: exclude {patterns} (case insensitive)"
