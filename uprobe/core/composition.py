"""Base-composition attributes: base fractions, homopolymer runs and windowed base counts.

They express sequence rules such as those of Xu et al. 2009 (PNAS 106:2289,
doi:10.1073/pnas.0812506106), which Oligostan applies to smiFISH probes (Tsanov et al. 2016,
NAR 44:e165, doi:10.1093/nar/gkw784): A < 28 %, no AAAA, C 22-28 %, no CCCC, and no 4 C in
any 6 consecutive bases within the first 12. Each is an attribute, filtered like any other; the
filters run before remove_overlap, so overlap removal chooses among passing windows.
"""
import math

COMPOSITION_TYPES = ('base_fraction', 'max_homopolymer', 'max_window_count')
# Options of each type; `bases` is required except for max_homopolymer, `window` for max_window_count.
COMPOSITION_OPTIONS = {
    'base_fraction': ('bases', 'positions'),
    'max_homopolymer': ('bases', 'positions'),
    'max_window_count': ('bases', 'window', 'positions'),
}
DNA_BASES = 'ACGT'


def _bases(value, kind, where):
    if value is None:
        if kind == 'max_homopolymer':
            return DNA_BASES
        raise ValueError(f'{where}.bases is required for {kind}, e.g. bases: A (or GC for G + C)')
    if isinstance(value, list) and all(isinstance(base, str) for base in value):
        value = ''.join(value)
    if not isinstance(value, str) or not value:
        raise ValueError(f'{where}.bases must be a string of bases such as A or GC, got {value!r}')
    bases = value.upper().replace('U', 'T')
    wrong = sorted(set(bases) - set(DNA_BASES))
    if wrong:
        raise ValueError(f'{where}.bases may contain only A, C, G, T (or U), got {value!r}')
    return ''.join(dict.fromkeys(bases))


def _positions(value, where):
    if value is None:
        return None
    hint = 'positions: [1, 12]   # 1-based and inclusive from the 5\' end; negative counts from the 3\' end (-1 = last)'
    if not (isinstance(value, list) and len(value) == 2
            and all(isinstance(p, int) and not isinstance(p, bool) and p != 0 for p in value)):
        raise ValueError(f'{where}.positions must be two non-zero integers [start, end], got {value!r}; {hint}')
    start, end = value
    if (start > 0) == (end > 0) and start > end:
        raise ValueError(f'{where}.positions: start {start} is after end {end}; {hint}')
    return start, end


def composition_settings(kind, attr, where='attribute'):
    """The checked options of a composition attribute: {'bases', 'positions', 'window'}; ValueError if invalid."""
    settings = {'bases': _bases(attr.get('bases'), kind, where), 'positions': _positions(attr.get('positions'), where)}
    if kind == 'max_window_count':
        window = attr.get('window')
        if not isinstance(window, int) or isinstance(window, bool) or window < 1:
            raise ValueError(f'{where}.window must be a positive integer (the window length in nt), got {window!r}')
        settings['window'] = window
    return settings


def _part(seq, positions):
    """The bases of ``seq`` at 1-based inclusive ``positions`` (negative: from the 3' end), clipped to the sequence."""
    seq = seq.upper().replace('U', 'T')
    if positions is None:
        return seq
    n = len(seq)
    start, end = (p if p > 0 else n + p + 1 for p in positions)
    return seq[max(start, 1) - 1:max(end, 0)]


def base_fraction(seq, bases, positions=None):
    """Fraction (0-1) of the bases of ``seq`` (within ``positions``) that are one of ``bases``."""
    if not isinstance(seq, str):
        return math.nan
    part = _part(seq, positions)
    if not part:
        return math.nan
    return sum(part.count(base) for base in bases) / len(part)


def max_homopolymer(seq, bases=DNA_BASES, positions=None):
    """Length of the longest run of one base of ``bases`` in ``seq`` (within ``positions``); 0 if none."""
    if not isinstance(seq, str):
        return math.nan
    longest = run = 0
    previous = None
    for base in _part(seq, positions):
        run = run + 1 if base == previous else 1
        previous = base
        if base in bases and run > longest:
            longest = run
    return longest


def max_window_count(seq, bases, window, positions=None):
    """Highest number of ``bases`` in any ``window`` consecutive bases of ``seq`` that lie within ``positions``.

    A part shorter than the window is counted whole.
    """
    if not isinstance(seq, str):
        return math.nan
    hits = [base in bases for base in _part(seq, positions)]
    if len(hits) <= window:
        return sum(hits)
    count = best = sum(hits[:window])
    for i in range(window, len(hits)):
        count += hits[i] - hits[i - window]
        best = max(best, count)
    return best


def composition_function(kind, settings):
    """A function of one sequence that computes the attribute ``kind`` with ``settings``."""
    bases, positions = settings['bases'], settings['positions']
    if kind == 'base_fraction':
        return lambda seq: base_fraction(seq, bases, positions)
    if kind == 'max_homopolymer':
        return lambda seq: max_homopolymer(seq, bases, positions)
    if kind == 'max_window_count':
        return lambda seq: max_window_count(seq, bases, settings['window'], positions)
    raise ValueError(f'not a composition attribute: {kind}')
