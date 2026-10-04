import json
import re
import pandas as pd
from uprobe.core.report.compact import build_report_data, save_compact_report


def fixture_data():
    raw = pd.DataFrame({"probe_id": ["p1", "p2", "p3"], "target": ["A", "A", "B"], "score": [0, 2, 9], "custom.part1": ["ACGT", "TGCA", "AAAA"]})
    protocol = {"name": "Custom", "targets": ["A", {"B": "AAAA"}], "attributes": {"score": {"type": "self_match"}}, "post_process": {"filters": {"score": {"condition": "score <= 4"}}}, "summary": {"report": {"key_metrics": ["score"], "table_columns": ["target", "score"]}}}
    return raw, protocol


def test_actual_counts_partial_coverage_and_custom_columns():
    raw, protocol = fixture_data()
    data = build_report_data(raw.iloc[:1], protocol, raw)
    assert (data["rawCount"], data["filterCount"], data["finalCount"]) == (3, 2, 1)
    assert data["covered"] == 1 and data["targetCount"] == 2
    assert data["state"] == "partial"
    assert data["tableColumns"] == list(raw.columns)
    assert data["statistics"][0]["median"] == 0
    assert data["diagnostics"][0]["rejected"] == 1


def test_empty_result_retains_raw_and_diagnostics(tmp_path):
    raw, protocol = fixture_data()
    data = build_report_data(raw.iloc[:0], protocol, raw)
    assert data["state"] == "empty" and data["finalCount"] == 0
    assert len(data["rawRows"]) == 3 and data["rows"] == []
    assert save_compact_report(raw.iloc[:0], protocol, tmp_path / "empty.html", raw)


def test_standalone_counts_are_unknown():
    raw, protocol = fixture_data()
    data = build_report_data(raw, protocol)
    assert data["rawCount"] is None and data["filterCount"] is None
    assert data["diagnostics"][0]["passed"] is None


def test_no_filters_not_reported_as_full_pass_rate():
    raw, protocol = fixture_data()
    protocol["post_process"] = {}
    data = build_report_data(raw, protocol, raw)
    assert data["filtersConfigured"] is False and data["filterCount"] is None


def test_safe_offline_english_report_and_exact_export_data(tmp_path):
    raw, protocol = fixture_data()
    protocol["name"] = '</script><script>alert("x")</script>'
    raw.loc[0, "custom.part1"] = '<img src=x onerror=alert(1)>'
    path = save_compact_report(raw, protocol, tmp_path / "safe.html", raw)
    html = path.read_text()
    assert '<html lang="en">' in html
    assert 'Run information' not in html
    assert 'View fields' not in html
    assert 'All fields are shown.' in html
    assert not re.search(r"[\u4e00-\u9fff]", html)
    assert 'src="https://' not in html and 'href="https://' not in html
    assert '<script>alert("x")</script>' not in html
    payload = re.search(r'<script type="application/json" id="report-data">(.*?)</script>', html, re.S).group(1)
    data = json.loads(payload)
    assert data["rows"][0]["score"] == 0
    assert data["rows"][0]["custom.part1"] == '<img src=x onerror=alert(1)>'
