"""Machine-readable run outputs for scripts and agents.

``summary.json`` states in one small file whether a run produced probes, how
many candidates each target kept after every step, and which step removed
them. The final probe table is also written as plain TSV, which is far easier
to parse than XLSX. Neither file contains per-probe data beyond the TSV itself.
"""
import csv
import json
import math
import os
import uuid
from pathlib import Path

import pandas as pd

from uprobe import __version__
from ..hybridization import effective_conditions
from ..offtarget import OFF_TARGETS_SUFFIX, parse_gene_hits
from ..process.equal_space import number_desired

SUMMARY_FILENAME = "summary.json"
SUMMARY_SCHEMA_VERSION = 1

# Early-exit statuses are set by the workflow; the others are derived here.
STATUS_MESSAGES = {
    "no_valid_targets": "No valid targets remained after validation against the genome annotation.",
    "no_target_sequences": "No target region sequences were extracted.",
    "no_probes_constructed": "No probes were constructed from the target regions.",
}
# A run with one of these statuses finished without any probe (`uprobe run` exits with 3).
NO_PROBE_STATUSES = ("no_probes", *STATUS_MESSAGES)


def _tsv_cell(value):
    if isinstance(value, set):
        value = sorted(value, key=str)
    if isinstance(value, (list, tuple, dict)):
        value = json.dumps(value, default=str)
    if isinstance(value, str):
        # Plain TSV has no quoting: tabs and line breaks inside a value become spaces.
        return value.replace("\t", " ").replace("\r", " ").replace("\n", " ")
    return value


def save_tsv(df: pd.DataFrame, path) -> Path:
    """Write every column as unquoted tab-separated text; list/dict cells are JSON."""
    table = df.copy()
    for index, dtype in enumerate(table.dtypes):
        if dtype == object or pd.api.types.is_string_dtype(dtype):
            table.iloc[:, index] = table.iloc[:, index].map(_tsv_cell)
    table.to_csv(path, sep="\t", index=False, lineterminator="\n", quoting=csv.QUOTE_NONE)
    return Path(path)


def _target_column(df):
    return next((c for c in ("target", "gene") if c in df.columns), None)


def _target_counts(df):
    """Rows per target, or None when the step was not reached."""
    if df is None:
        return None
    column = _target_column(df)
    if column is None:
        return {}
    return {str(k): int(v) for k, v in df[column].dropna().astype(str).value_counts(sort=False).items()}


def _built_rows(df, probe_names):
    """Rows whose top-level probe sequences were all built (failed nodes are '')."""
    if df is None:
        return None
    columns = [name for name in probe_names if name in df.columns]
    if not columns:
        return df
    built = pd.Series(True, index=df.index)
    for column in columns:
        built &= df[column].fillna("").astype(str).str.len() > 0
    return df[built]


def _number(value):
    value = float(value)
    return round(value, 4) if math.isfinite(value) else None


def _attribute_stats(df, protocol):
    """min/median/max of summary.attributes (or all numeric attributes) on the final table."""
    if df is None:
        return {}
    configured = (protocol.get("summary") or {}).get("attributes") or []
    names = [name for name in (configured or (protocol.get("attributes") or {})) if isinstance(name, str)]
    stats = {}
    for name in names:
        if name not in df.columns or isinstance(df[name], pd.DataFrame):
            continue
        series = df[name]
        if not pd.api.types.is_numeric_dtype(series) or pd.api.types.is_bool_dtype(series):
            continue
        values = pd.to_numeric(series, errors="coerce").dropna()
        values = values[values.map(math.isfinite)]
        stats[name] = {
            "count": int(len(values)),
            "min": _number(values.min()) if len(values) else None,
            "median": _number(values.median()) if len(values) else None,
            "max": _number(values.max()) if len(values) else None,
        }
    return stats


def _hybridization_conditions(protocol):
    """Tm/fold_score conditions in effect per attribute, or None if there are none or they are invalid."""
    try:
        return effective_conditions(protocol) or None
    except Exception:
        return None


def _requested_probes(protocol, target):
    """post_process.equal_space.number_desired for one target, or None."""
    post_process = protocol.get("post_process")
    value = number_desired(post_process.get("equal_space"), target) if isinstance(post_process, dict) else None
    return value if type(value) is int and value > 0 else None


def _off_target_summary(report, final):
    """The expression table in use, per mapped_genes attribute the off-target genes found, and the pool limit.

    ``report`` is what add_attributes recorded; final_* counts come from the
    <attribute>_off_targets column of the final table, and ``pool`` from
    post_process.off_target_pool.
    """
    if not isinstance(report, dict) or not report.get("attributes"):
        return None
    attributes = {}
    for name, stats in report["attributes"].items():
        entry = {key: value for key, value in stats.items() if key != "missing_labels"}
        column = name + OFF_TARGETS_SUFFIX
        if final is not None and column in final.columns:
            genes = {label for cell in final[column] for label in parse_gene_hits(cell)}
            entry["final_off_target_genes"] = len(genes)
            if "missing_labels" in stats:
                entry["final_missing_from_table"] = sorted(genes & set(stats["missing_labels"]))
        attributes[name] = entry
    pool = final.attrs.get("off_target_pool") if final is not None else None
    return {"gene_expression": report.get("gene_expression"), "attributes": attributes, "pool": pool}


def _relative(path, output_dir):
    if not path:
        return None
    path = Path(path)
    if output_dir is not None:
        if not path.is_absolute():
            path = Path(output_dir) / path
        if not path.exists():
            return None
        try:
            return path.relative_to(output_dir).as_posix()
        except ValueError:
            return str(path)
    return path.name


def build_run_summary(protocol: dict, run: dict, output_dir=None) -> dict:
    """Assemble the summary from facts the workflow recorded in ``run``.

    ``run`` keys (all optional): started, finished (datetimes), step, status,
    error, requested_targets, target_validation ({'valid', 'invalid'}),
    targets (extracted windows), combined (targets + constructed probes),
    raw (after attributes), final (after post-processing), files
    ({'xlsx', 'raw_xlsx', 'tsv', 'html'}), no_filtered_probes, raw_file,
    genome_source (the genome entry's ``source``, as written by ``uprobe genome fetch``),
    off_targets (the off-target statistics add_attributes recorded).
    """
    protocol = protocol if isinstance(protocol, dict) else {}
    region = (protocol.get("extracts") or {}).get("target_region") or {}
    requested = [str(t) for t in run.get("requested_targets") or []]
    validation = run.get("target_validation")
    if validation:
        valid, invalid = [str(t) for t in validation["valid"]], [str(t) for t in validation["invalid"]]
    elif run.get("targets") is not None:
        # Validation is skipped for genome regions: every requested target was used.
        valid, invalid = list(requested), []
    else:
        # The run stopped before targets were validated.
        valid, invalid = None, []

    final = run.get("final")
    stages = list(final.attrs.get("post_process_stages", [])) if final is not None else []
    stage_counts = list(final.attrs.get("post_process_target_counts", [])) if final is not None else []
    if len(stage_counts) != len(stages):
        stage_counts = []
    windows = _target_counts(run.get("targets"))
    built_rows = _built_rows(run.get("combined"), list(protocol.get("probes") or {}))
    constructed = _target_counts(built_rows)
    before_post = _target_counts(run.get("raw"))
    final_counts = _target_counts(final)

    names = list(dict.fromkeys((requested if valid is None else valid) + list(windows or {})))
    targets = {}
    eliminated = {}
    for name in names:
        entry = {
            "candidate_windows": windows.get(name, 0) if windows is not None else None,
            "probes_constructed": constructed.get(name, 0) if constructed is not None else None,
            "after_stage": {},
            "final": final_counts.get(name, 0) if final_counts is not None else 0,
            "requested_probes": _requested_probes(protocol, name),
            "top_removal_stage": None,
        }
        if before_post is not None and stage_counts:
            before, most = before_post.get(name, 0), 0
            for stage, counts in zip(stages, stage_counts):
                after = counts.get(name, 0)
                entry["after_stage"][stage["stage"]] = after
                if before - after > most:
                    most, entry["top_removal_stage"] = before - after, stage["stage"]
                if before > 0 and after == 0:
                    eliminated[name] = stage["stage"]
                before = after
        targets[name] = entry

    error = run.get("error")
    final_total = len(final) if final is not None else 0
    without_probes = [name for name, entry in targets.items() if entry["final"] == 0]
    # Targets that kept fewer probes than equal_space requested (all passing ones are kept).
    shortfall = {name: {"requested": entry["requested_probes"], "final": entry["final"]}
                 for name, entry in targets.items()
                 if final is not None and entry["requested_probes"] and entry["final"] < entry["requested_probes"]}
    if error is not None:
        status = "error"
    elif run.get("status"):
        status = run["status"]
    elif final_total == 0:
        status = "no_probes"
    elif invalid or without_probes or shortfall:
        status = "partial"
    else:
        status = "ok"

    files = run.get("files") or {}
    output_files = {key: _relative(files.get(key), output_dir) for key in ("xlsx", "raw_xlsx", "tsv", "html")}
    output_files["summary"] = SUMMARY_FILENAME
    raw_file = run.get("raw_file")

    warnings = []
    if invalid:
        warnings.append({"code": "invalid_targets", "targets": invalid,
                         "message": f"{len(invalid)} target(s) not found in the genome annotation and skipped."})
    if windows is not None:
        no_windows = [name for name in names if not windows.get(name)]
        if no_windows:
            warnings.append({"code": "no_candidate_windows", "targets": no_windows,
                             "message": f"{len(no_windows)} target(s) yielded no candidate windows "
                                        "(sequence missing or shorter than extracts.target_region.length)."})
    if built_rows is not None and len(built_rows) < len(run["combined"]):
        failed = len(run["combined"]) - len(built_rows)
        warnings.append({"code": "probe_construction_failed", "count": failed,
                         "message": f"{failed} candidate(s) have an empty probe sequence; check the probe expressions."})
    for stage in stages:
        if stage.get("before", 0) > 0 and stage.get("after") == 0:
            warnings.append({"code": "stage_removed_all", "stage": stage["stage"], "removed": stage["before"],
                             "message": f"Post-process stage '{stage['stage']}' removed all {stage['before']} remaining candidates."})
    if final is not None and without_probes:
        warnings.append({"code": "targets_without_probes", "targets": without_probes,
                         "eliminated_by": {name: eliminated.get(name) for name in without_probes},
                         "message": f"{len(without_probes)} target(s) have no final probes."})
    if shortfall:
        warnings.append({"code": "fewer_probes_than_requested", "targets": list(shortfall), "counts": shortfall,
                         "message": f"{len(shortfall)} target(s) have fewer probes than requested "
                                    "(post_process.equal_space.number_desired); all their passing probes were kept."})
    off_targets = _off_target_summary(run.get("off_targets"), final)
    missing = {name: entry["missing_from_table"] for name, entry in (off_targets or {}).get("attributes", {}).items()
               if entry.get("missing_from_table")}
    if missing:
        warnings.append({"code": "off_target_genes_not_in_expression_table", "attributes": missing,
                         "message": f"{sum(missing.values())} off-target gene(s) are not in the gene_expression "
                                    "table and were counted as expressed (see off_targets)."})

    if status == "error":
        message = f"Run failed during {run.get('step')}: {error}"
    elif status in STATUS_MESSAGES:
        message = STATUS_MESSAGES[status]
    elif status == "no_probes":
        message = "No candidate passed post-processing." + (f" Raw candidates: {raw_file}." if raw_file else "")
    elif status == "partial":
        gaps = []
        if without_probes or invalid:
            gaps.append(f"{len(without_probes) + len(invalid)} of {len(requested) or len(names)} "
                        "requested target(s) have no probes")
        fewer = [name for name, counts in shortfall.items() if counts["final"] > 0]
        if fewer:
            gaps.append(f"{len(fewer)} target(s) have fewer probes than requested")
        message = f"{final_total} probes designed; " + "; ".join(gaps) + "."
    else:
        message = f"{final_total} probes designed for all {len(names)} target(s)."

    started, finished = run.get("started"), run.get("finished")
    genome_source = run.get("genome_source")
    return {
        "schema_version": SUMMARY_SCHEMA_VERSION,
        "uprobe_version": __version__,
        "protocol_name": protocol.get("name"),
        "status": status,
        "message": message,
        "stopped_at": run.get("step") if status in STATUS_MESSAGES or status == "error" else None,
        "error": f"{type(error).__name__}: {error}" if error is not None else None,
        "started_at": started.isoformat(timespec="seconds") if started else None,
        "finished_at": finished.isoformat(timespec="seconds") if finished else None,
        "duration_seconds": round((finished - started).total_seconds(), 3) if started and finished else None,
        "inputs": {
            "genome": protocol.get("genome"),
            "extract_source": region.get("source"),
            "target_length": region.get("length"),
            "step": region.get("step"),
            "requested_targets": requested,
            "valid_targets": valid,
            "invalid_targets": invalid,
            "hybridization": protocol.get("hybridization"),
            "hybridization_conditions": _hybridization_conditions(protocol),
            "genome_source": genome_source if isinstance(genome_source, dict) else None,
        },
        "counts": {
            "requested_targets": len(requested),
            "valid_targets": len(valid) if valid is not None else None,
            "candidate_windows": len(run["targets"]) if run.get("targets") is not None else None,
            "probes_constructed": len(built_rows) if built_rows is not None else None,
            "final_probes": final_total,
            "targets_with_probes": sum(entry["final"] > 0 for entry in targets.values()),
        },
        "stages": stages,
        "targets": targets,
        "attributes": _attribute_stats(final, protocol),
        "off_targets": off_targets,
        "output_files": output_files,
        "flags": {"no_filtered_probes": bool(run.get("no_filtered_probes")), "raw_file": raw_file},
        "warnings": warnings,
    }


def summary_json(summary: dict) -> str:
    return json.dumps(summary, indent=2, ensure_ascii=False, default=str, allow_nan=False)


def write_run_summary(summary: dict, path) -> Path:
    """Write atomically so a polling reader never sees a partial file."""
    path = Path(path)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        temporary.write_text(summary_json(summary) + "\n", encoding="utf-8")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    return path
