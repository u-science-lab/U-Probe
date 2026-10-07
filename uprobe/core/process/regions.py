import pandas as pd


def region_bounds(df: pd.DataFrame) -> pd.DataFrame:
    """Return numeric 1-based `start`/`end` for each candidate, aligned to df.index.

    Genome windows only carry `sub_region` ("101-140"), exon/UTR windows write
    "101_140", and direct sequences carry both forms; explicit columns win.
    """
    bounds = pd.DataFrame(index=df.index, columns=['start', 'end'], dtype='float')
    if 'sub_region' in df.columns:
        parsed = df['sub_region'].astype(str).str.extract(r'^(\d+)[-_](\d+)$')
        bounds['start'] = pd.to_numeric(parsed[0])
        bounds['end'] = pd.to_numeric(parsed[1])
    for column in ('start', 'end'):
        if column in df.columns:
            explicit = pd.to_numeric(df[column], errors='coerce')
            bounds[column] = explicit.fillna(bounds[column])
    if bounds.isna().any().any():
        bad = df.index[bounds.isna().any(axis=1)][:5].tolist()
        raise ValueError(f"Cannot determine start/end for candidate rows {bad}")
    return bounds.astype(int)
