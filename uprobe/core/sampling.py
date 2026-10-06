"""Shared validation and window generation for all target sources."""


def normalize_sampling(length, step=None, overlap=None):
    if type(length) is int:
        lower = upper = length
    elif isinstance(length, (list, tuple)) and len(length) == 2:
        lower, upper = length
    else:
        raise ValueError('length must be a positive integer or [min, max]')
    if any(type(v) is not int or v <= 0 for v in (lower, upper)) or lower > upper:
        raise ValueError('length bounds must be positive integers with min <= max')
    if step is None:
        if overlap is not None:
            if type(length) is not int:
                raise ValueError('range length requires step; legacy overlap supports fixed length only')
            if type(overlap) is not int or not 0 <= overlap < lower:
                raise ValueError('overlap must be an integer with 0 <= overlap < length')
            step = lower - overlap
        else:
            step = 1
    if type(step) is not int or step <= 0:
        raise ValueError('step must be a positive integer')
    return lower, upper, step


def iter_target_windows(seq, length, step=None, overlap=None):
    """Yield zero-based start, exclusive end and sequence; ranges are inclusive."""
    lower, upper, step = normalize_sampling(length, step, overlap)
    for start in range(0, len(seq) - lower + 1, step):
        for size in range(lower, min(upper, len(seq) - start) + 1):
            yield start, start + size, seq[start:start + size]
