import pandas as pd
import pytest
from zipfile import ZipFile
from uprobe.core import api as module


@pytest.mark.parametrize("raw_csv", [True, False])
def test_empty_filters_preserve_raw_without_returning_it(tmp_path, monkeypatch, raw_csv):
    api = object.__new__(module.UProbeAPI)
    api.protocol = {"name": "empty_test", "post_process": {"filters": {"score": {"condition": "score < 0"}}}}
    api.output_dir = tmp_path
    api.genome = {}
    raw = pd.DataFrame({"probe_id": ["p1"], "score": [10]})
    monkeypatch.setattr(module, "add_attributes", lambda *args: raw.copy())
    monkeypatch.setattr(module, "post_process", lambda *args: pd.DataFrame())
    result = api.post_process_probes(raw, raw_csv=raw_csv)
    assert result.empty
    assert api.no_filtered_probes is True
    assert api.raw_file.endswith(".xlsx")
    assert not list(tmp_path.glob("*.csv"))
    with ZipFile(tmp_path / api.raw_file) as book:
        assert b"p1" in book.read("xl/sharedStrings.xml")
    with ZipFile(tmp_path / api._csv_filename) as book:
        assert b'r="2"' not in book.read("xl/worksheets/sheet1.xml")


def test_passing_probes_clear_previous_warning(tmp_path, monkeypatch):
    api = object.__new__(module.UProbeAPI)
    api.protocol = {"name": "pass_test", "post_process": {"filters": {"score": {"condition": "score < 20"}}}}
    api.output_dir = tmp_path
    api.genome = {}
    api.no_filtered_probes = True
    api.raw_file = "old_raw.csv"
    raw = pd.DataFrame({"probe_id": ["p1"], "score": [10]})
    monkeypatch.setattr(module, "add_attributes", lambda *args: raw.copy())
    monkeypatch.setattr(module, "post_process", lambda *args: raw.copy())
    assert len(api.post_process_probes(raw)) == 1
    assert api.no_filtered_probes is False
    assert api.raw_file is None
