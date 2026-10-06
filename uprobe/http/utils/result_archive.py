"""Keep downloadable task archives limited to reports and probe tables."""

from pathlib import Path
import os
import shutil
import tempfile
import zipfile
import xml.etree.ElementTree as ET

RESULT_SUFFIXES = frozenset({".xlsx", ".html"})


def _has_table_rows(source, entry):
    """Inspect worksheet data, preserving unreadable legacy files conservatively."""
    with tempfile.TemporaryFile() as stream:
        with source.open(entry) as reader:
            shutil.copyfileobj(reader, stream)
        stream.seek(0)
        try:
            with zipfile.ZipFile(stream) as workbook:
                for name in workbook.namelist():
                    if not name.startswith('xl/worksheets/') or not name.endswith('.xml'):
                        continue
                    with workbook.open(name) as sheet:
                        for _, row in ET.iterparse(sheet, events=('end',)):
                            if row.tag != '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}row':
                                continue
                            if int(row.get('r', '1')) > 1 and any(
                                node.text for node in row.iter()
                                if node.tag.rsplit('}', 1)[-1] in ('v', 't', 'f')
                            ):
                                return True
                            row.clear()
                return False
        except (zipfile.BadZipFile, ET.ParseError, ValueError):
            return True


def restrict_result_archive(archive_path: Path) -> None:
    """Remove configuration/log entries and empty tables, atomically replacing the ZIP."""
    temporary_path = None
    try:
        with zipfile.ZipFile(archive_path) as source:
            entries = source.infolist()
            allowed = [entry for entry in entries
                       if not entry.is_dir()
                       and Path(entry.filename).suffix in RESULT_SUFFIXES
                       and (Path(entry.filename).suffix != '.xlsx'
                            or _has_table_rows(source, entry))]
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
