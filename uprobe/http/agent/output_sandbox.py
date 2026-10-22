"""Shared path contract for Web and CLI agent sessions."""

from __future__ import annotations

import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal

from uprobe.http.utils.agent_store import safe_path_component
from uprobe.http.utils.paths import get_server_root, get_data_dir, get_probes_dir, get_genomes_yaml


AgentEntrypoint = Literal["cli", "web"]
_DEFAULT_CLI_SESSION_ID: str | None = None


@dataclass(frozen=True)
class AgentRuntimePaths:
    """Canonical directories and resource files exposed to agent runtimes."""

    entrypoint: AgentEntrypoint
    workspace_root: Path
    project_root: Path
    output_root: Path
    memory_root: Path
    data_dir: Path
    probes_dir: Path
    genomes_path: Path


def cli_user_key(user: str | None = None) -> str:
    return safe_path_component(
        (user or "").strip()
        or (os.environ.get("UPROBE_AGENT_USER") or "").strip()
        or os.environ.get("USER", "")
        or os.environ.get("USERNAME", "")
        or "cli",
    )


def _new_cli_session_id() -> str:
    return f"cli_{datetime.now(timezone.utc).strftime('%Y%m%d%H%M%S%f')}"


def cli_session_component(chat_id: str | None = None) -> str:
    global _DEFAULT_CLI_SESSION_ID

    manual = (os.environ.get("UPROBE_AGENT_SESSION") or "").strip()
    if manual:
        s = safe_path_component(manual)
        if s:
            return s
    if chat_id and str(chat_id).strip():
        s = safe_path_component(str(chat_id).strip())
        if s:
            return s
    if _DEFAULT_CLI_SESSION_ID is None:
        _DEFAULT_CLI_SESSION_ID = _new_cli_session_id()
    return _DEFAULT_CLI_SESSION_ID


def _resolve_project_root(workspace: Path) -> Path:
    ws = workspace.expanduser().resolve()
    if (ws / "data" / "genomes.yaml").exists():
        return ws
    return get_server_root().expanduser().resolve()


def _resolve_output_root(
    *,
    entrypoint: AgentEntrypoint,
    workspace: Path,
    output_dir: Path | None,
    user: str | None,
    session_id: str | None,
) -> Path:
    if output_dir is not None:
        return output_dir.expanduser().resolve()
    preset = os.environ.get("UPROBE_OUTPUT_DIR", "").strip()
    if preset:
        return Path(preset).expanduser().resolve()
    if entrypoint == "web":
        raise ValueError("Web agent sessions must provide output_dir")

    return (
        workspace.expanduser().resolve()
        / "outputs"
        / "agent"
        / cli_user_key(user)
        / cli_session_component(session_id)
    ).resolve()


def _resolve_memory_root(
    *,
    entrypoint: AgentEntrypoint,
    workspace: Path,
    memory_dir: Path | None,
    session_id: str | None,
) -> Path:
    if memory_dir is not None:
        return memory_dir.expanduser().resolve()
    if entrypoint == "web":
        raise ValueError("Web agent sessions must provide memory_dir")
    return (
        workspace.expanduser().resolve()
        / ".pantheon"
        / "memory"
        / cli_session_component(session_id)
    ).resolve()


def resolve_agent_paths(
    *,
    entrypoint: AgentEntrypoint,
    workspace: Path | None = None,
    output_dir: Path | None = None,
    memory_dir: Path | None = None,
    user: str | None = None,
    session_id: str | None = None,
) -> AgentRuntimePaths:
    """Resolve the only path contract agent code should consume."""

    workspace_root = (workspace or get_server_root()).expanduser().resolve()
    project_root = _resolve_project_root(workspace_root)
    data_dir = get_data_dir()
    effective_session_id = cli_session_component(session_id) if entrypoint == "cli" else session_id
    output_root = _resolve_output_root(
        entrypoint=entrypoint,
        workspace=workspace_root,
        output_dir=output_dir,
        user=user,
        session_id=effective_session_id,
    )
    memory_root = _resolve_memory_root(
        entrypoint=entrypoint,
        workspace=workspace_root,
        memory_dir=memory_dir,
        session_id=effective_session_id,
    )

    output_root.mkdir(parents=True, exist_ok=True)
    memory_root.mkdir(parents=True, exist_ok=True)

    return AgentRuntimePaths(
        entrypoint=entrypoint,
        workspace_root=workspace_root,
        project_root=project_root,
        output_root=output_root,
        memory_root=memory_root,
        data_dir=data_dir,
        probes_dir=get_probes_dir(),
        genomes_path=get_genomes_yaml(),
    )


def apply_agent_runtime_env(paths: AgentRuntimePaths) -> None:
    """Expose canonical project/resource/output paths to agent subprocesses."""

    os.environ["UPROBE_PROJECT_ROOT"] = str(paths.project_root)
    os.environ["UPROBE_WORKSPACE_ROOT"] = str(paths.workspace_root)
    os.environ["UPROBE_DATA_DIR"] = str(paths.data_dir)
    os.environ.pop("UPROBE_PROBE_JSON", None)
    os.environ["UPROBE_PROBES_DIR"] = str(get_probes_dir())
    os.environ["UPROBE_GENOMES_PATH"] = str(paths.genomes_path)
    os.environ["UPROBE_OUTPUT_DIR"] = str(paths.output_root)


def probe_resource_context() -> str:
    """Update old conversations without reading templates on a greeting."""
    return (
        "[Application probe resource context]\n"
        f"Current built-in library: {get_probes_dir()} (<method>.yaml). "
        "This replaces earlier probe.json resource instructions. For method lookups "
        "use Python uprobe.http.utils.probe_library.load_builtin_probes(); it returns "
        "method names mapped to configurations, exactly like the web workflow API. "
        "Immediately before a design, call load_builtin_probe(method_name) from the "
        "same module to read the latest configuration. Never use remembered presets "
        "or the retired JSON file. A file-manager-only agent delegates this lookup "
        "to Probe_Designer. A greeting does not require a lookup.\n"
        "[/Application probe resource context]"
    )


def resolve_cli_agent_output_dir(chat_id: str | None) -> Path:
    return resolve_agent_paths(entrypoint="cli", session_id=chat_id).output_root


def apply_cli_agent_runtime_env(*, workspace: Path, chat_id: str | None) -> Path:
    """
    Backward-compatible CLI helper.

    New code should call ``resolve_agent_paths`` and ``apply_agent_runtime_env`` directly.
    """

    paths = resolve_agent_paths(entrypoint="cli", workspace=workspace, session_id=chat_id)
    apply_agent_runtime_env(paths)
    return paths.output_root
