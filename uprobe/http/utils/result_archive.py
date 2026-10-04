"""Keep downloadable task archives limited to reports and probe tables."""

from pathlib import Path
import os
import shutil
import tempfile
import zipfile

RESULT_SUFFIXES = frozenset({".csv", ".html"})


def restrict_result_archive(archive_path: Path) -> None:
    """Remove legacy configuration/log entries, atomically replacing the ZIP."""
    temporary_path = None
    try:
        with zipfile.ZipFile(archive_path) as source:
            entries = source.infolist()
            allowed = [entry for entry in entries
                       if not entry.is_dir()
                       and Path(entry.filename).suffix in RESULT_SUFFIXES]
            if len(allowed) == len(entries):
                return
            with tempfile.NamedTemporaryFile(
                dir=archive_path.parent, suffix=".zip", delete=False
            ) as temporary:
                temporary_path = Path(temporary.name)
            with zipfile.ZipFile(temporary_path, "w", zipfile.ZIP_DEFLATED) as target:
                for entry in allowed:
                    with source.open(entry) as reader, target.open(entry, "w") as writer:
                        shutil.copyfileobj(reader, writer)
        os.replace(temporary_path, archive_path)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
