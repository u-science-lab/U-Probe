"""Presentation-only normalization; internal coordinates remain available."""
import ast
import pandas as pd


def format_transcripts(value):
    if isinstance(value, str):
        stripped = value.strip()
        if stripped.startswith(("[", "(")):
            try:
                value = ast.literal_eval(stripped)
            except (ValueError, SyntaxError):
                return value
    if isinstance(value, (list, tuple, set)):
        return "; ".join(str(item) for item in (sorted(value) if isinstance(value, set) else value))
    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return ""
    return str(value)


def prepare_result_table(df):
    result = df.drop(columns=["start", "end"], errors="ignore").copy()
    if "transcript_names" in result:
        result["transcript_names"] = result["transcript_names"].map(format_transcripts)
    return result
