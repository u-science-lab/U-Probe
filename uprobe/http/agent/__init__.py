"""
Chat agent of the U-Probe web server: an optional component.

The ``/agent`` endpoints (``uprobe.http.routers.agent``) run a PantheonOS agent
team (``session_manager``, ``output_sandbox``, ``templates/``). PantheonOS is not
a dependency of ``pip install uprobe``; the ``web-agent`` extra installs it::

    pip install 'uprobe[web-agent]'

Without it the server starts normally and the ``/agent`` endpoints answer 503.
This ``__init__`` imports neither PantheonOS nor the submodules, and nothing under
``uprobe.core`` imports this package.
"""

from __future__ import annotations

import importlib

WEB_AGENT_EXTRA = "web-agent"
WEB_AGENT_INSTALL = f"pip install 'uprobe[{WEB_AGENT_EXTRA}]'"
WEB_AGENT_UNAVAILABLE = (
    "The web chat agent is not installed on this server. "
    f"Install the '{WEB_AGENT_EXTRA}' extra ({WEB_AGENT_INSTALL}) and restart the server."
)


def web_agent_available() -> bool:
    """Whether the agent runtime of the ``web-agent`` extra (PantheonOS) can be imported."""

    try:
        importlib.import_module("pantheon")
    except ImportError:
        return False
    return True


def agent_unavailable_router():
    """Stand-in for the ``/agent`` router when the extra is missing: every endpoint answers 503."""

    from fastapi import APIRouter, HTTPException

    router = APIRouter(prefix="/agent", tags=["agent"])

    @router.api_route(
        "/{path:path}",
        methods=["GET", "POST", "PUT", "PATCH", "DELETE"],
        include_in_schema=False,
    )
    async def web_agent_not_installed(path: str):
        raise HTTPException(status_code=503, detail=WEB_AGENT_UNAVAILABLE)

    return router
