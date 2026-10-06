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


def sequence_groups(protocol):
    """Use configured sequence names, preserving the protocol's order."""
    return list(dict.fromkeys(["target_region", *(["target_parts"] if protocol.get("extracts", {}).get("target_region", {}).get("layout") else []), *protocol.get("extracts", {}), *protocol.get("probes", {})]))


def column_group(column, protocol):
    attributes = protocol.get("attributes", {}) or {}
    config = attributes.get(column, {})
    target = config.get("target", column) if isinstance(config, dict) else column
    if not isinstance(target, str):
        return None
    return next((group for group in sorted(sequence_groups(protocol), key=len, reverse=True)
                 if target == group or target.startswith(group + ".")), None)


def ordered_columns(columns, protocol):
    """Keep metadata first, then each sequence, its parts, and its attributes."""
    columns = list(columns)
    groups = {column: column_group(column, protocol) for column in columns}
    ordered = [column for column in columns if groups[column] is None]
    attributes = protocol.get("attributes", {}) or {}
    for group in sequence_groups(protocol):
        members = [column for column in columns if groups[column] == group]
        ordered.extend(column for column in members if column == group)
        ordered.extend(column for column in members if column != group and column not in attributes)
        ordered.extend(column for column in members if column != group and column in attributes)
    return ordered


def prepare_result_table(df, protocol=None):
    result = df.drop(columns=["start", "end"], errors="ignore").copy()
    if "transcript_names" in result:
        result["transcript_names"] = result["transcript_names"].map(format_transcripts)
    if protocol is not None:
        result = result.loc[:, ordered_columns(result.columns, protocol)]
    return result
