from io import BytesIO
import xml.etree.ElementTree as ET
from zipfile import ZipFile
import pandas as pd
from uprobe.core.report.formatting import prepare_result_table
from uprobe.core.report.excel import workbook_bytes, column_styles
from uprobe.core.report.formatting import ordered_columns


def data():
    df = pd.DataFrame({"probe_id": ["p1"], "transcript_names": [["TX1", "TX2"]], "target_region": ["ACGT"], "pad_probe": ["CCTT"], "pad_probe.part1": ["CC"], "pad_tm": [0], "amp_probe": ["GGAA"], "amp_tm": [38.5], "start": [1], "end": [40]})
    protocol = {"probes": {"pad_probe": {}, "amp_probe": {}}, "attributes": {"pad_tm": {"target": "pad_probe.part1"}, "amp_tm": {"target": "amp_probe"}}}
    return df, protocol


def test_presentation_keeps_internal_coordinates():
    df, _ = data()
    formatted = prepare_result_table(df)
    assert "start" not in formatted and "end" not in formatted
    assert formatted.transcript_names.tolist() == ["TX1; TX2"]
    assert df.transcript_names.tolist() == [["TX1", "TX2"]]
    assert "start" in df and "end" in df
    df.transcript_names = ["['TX1', 'TX2']"]
    assert prepare_result_table(df).transcript_names.tolist() == ["TX1; TX2"]


def test_sequence_components_and_attributes_share_color_group():
    df, protocol = data()
    styles = column_styles(df.columns, protocol)
    assert {styles[c]["group"] for c in ["pad_probe", "pad_probe.part1", "pad_tm"]} == {"pad_probe"}
    assert {styles[c]["group"] for c in ["amp_probe", "amp_tm"]} == {"amp_probe"}
    assert styles["pad_tm"]["background"] == styles["pad_probe.part1"]["background"]
    assert styles["pad_probe"]["header"] == styles["pad_tm"]["header"]
    assert styles["pad_probe"]["header"] != styles["amp_probe"]["header"]


def test_generic_order_uses_configured_targets_without_losing_columns():
    protocol = {"extracts": {"capture": {}}, "probes": {"custom": {}, "custom.long": {}},
                "attributes": {"unrelated_name": {"target": "custom.part1"},
                               "capture_quality": {"target": "capture"},
                               "long_quality": {"target": "custom.long.part1"},
                               "unknown": {"target": "missing"}}}
    columns = ["long_quality", "unrelated_name", "custom.part1", "target", "custom.long.part1",
               "capture_quality", "custom.long", "unknown", "capture", "custom", "extra"]
    expected = ["target", "unknown", "extra", "capture", "capture_quality", "custom",
                "custom.part1", "unrelated_name", "custom.long", "custom.long.part1", "long_quality"]
    assert ordered_columns(columns, protocol) == expected
    df = pd.DataFrame([{column: column for column in columns}])
    assert prepare_result_table(df, protocol).columns.tolist() == expected
    assert df.columns.tolist() == columns
    assert ordered_columns(columns, {}) == columns
    styles = column_styles(columns, protocol)
    assert styles["long_quality"]["group"] == "custom.long"
    assert styles["unrelated_name"]["group"] == "custom"


def test_xlsx_values_styles_and_navigation():
    df, protocol = data()
    with ZipFile(BytesIO(workbook_bytes(df, protocol))) as z:
        ns = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
        strings = ET.fromstring(z.read("xl/sharedStrings.xml"))
        values = ["".join(el.itertext()) for el in strings]
        assert "start" not in values and "end" not in values
        assert "TX1; TX2" in values
        sheet = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
        assert sheet.find("s:autoFilter", ns) is not None
        assert sheet.find("s:sheetViews/s:sheetView/s:pane", ns) is not None
        cells = {cell.attrib["r"]: cell for cell in sheet.findall(".//s:c", ns)}
        assert cells["F2"].find("s:v", ns).text == "0"
        formats = ET.fromstring(z.read("xl/styles.xml")).find("s:cellXfs", ns)
        fill = lambda cell: formats[int(cells[cell].attrib["s"])].attrib["fillId"]
        assert fill("E2") == fill("F2")
        assert fill("F2") != fill("H2")

        borders = ET.fromstring(z.read("xl/styles.xml")).find("s:borders", ns)
        for cell in cells.values():
            border_id = int(formats[int(cell.attrib["s"])].attrib["borderId"])
            assert all(borders[border_id].find("s:" + edge, ns).attrib.get("style") == "thin" for edge in ("left", "right", "top", "bottom"))


# -- Excel limits: never fail the run or lose data -----------------------------------------------

def excel_api(tmp_path, monkeypatch, table, report_name=''):
    """A UProbeAPI whose attributes step returns ``table``; post-processing keeps every row."""
    from uprobe.core import api as module
    api = object.__new__(module.UProbeAPI)
    api.protocol = {"name": "big", "probes": {}, "post_process": {"sorts": {"is_ascending": ["score"]}},
                    "summary": {"report_name": report_name}}
    api.output_dir = tmp_path
    api.genome = {}
    api._generate_html, api._include_plots, api._embed_plots = True, False, False
    monkeypatch.setattr(module, "add_attributes", lambda *args: table.copy())
    return api


def read_cell(path, cell):
    import openpyxl
    return openpyxl.load_workbook(path, read_only=True)["Probes"][cell].value


def test_long_cells_are_cut_in_xlsx_and_kept_in_tsv(tmp_path, monkeypatch):
    from uprobe.core.report.excel import EXCEL_MAX_CELL_CHARS
    from uprobe.core.report.run_summary import build_run_summary
    hits = ", ".join(f"GENE{i}:2" for i in range(6000))  # a repeat window's off-target list
    assert len(hits) > EXCEL_MAX_CELL_CHARS
    table = pd.DataFrame({"probe_id": ["p1", "p2"], "target": ["G", "G"], "score": [1, 2],
                          "hits_off_targets": [hits, "ACTG1:0"]})
    api = excel_api(tmp_path, monkeypatch, table, report_name="rna_report")
    result = api.post_process_probes(table, raw_csv=True)
    files = api.output_files
    assert len(result) == 2 and files["raw_xlsx"] and files["xlsx"] and files["tsv"] and files["raw_tsv"]
    for xlsx, tsv in ((files["raw_xlsx"], files["raw_tsv"]), (files["xlsx"], files["tsv"])):
        value = read_cell(tmp_path / xlsx, "D2")
        assert len(value) == EXCEL_MAX_CELL_CHARS and value.startswith(hits[:1000])
        assert value.endswith(f"...[truncated: {len(hits):,} characters; full value in {tsv}]")
        full = pd.read_csv(tmp_path / tsv, sep="\t")
        assert full.loc[full.probe_id == "p1", "hits_off_targets"].item() == hits
    assert [entry["truncated_cells"] for entry in api.excel_limits] == [{"hits_off_targets": 1}] * 2
    # The HTML report embeds the tables as workbooks too.
    assert api.generate_report(result)["html_reports"]
    summary = build_run_summary({}, {"final": result, "files": files, "excel_limits": api.excel_limits}, tmp_path)
    warning = next(w for w in summary["warnings"] if w["code"] == "excel_limits")
    assert warning["tables"][0]["tsv"] == files["raw_tsv"] and summary["output_files"]["raw_tsv"] == files["raw_tsv"]


def test_tables_longer_than_a_worksheet_fall_back_to_tsv(tmp_path, monkeypatch):
    from uprobe.core.report import excel
    monkeypatch.setattr(excel, "EXCEL_MAX_ROWS", 4)  # a header and 3 data rows
    table = pd.DataFrame({"probe_id": [f"p{i}" for i in range(5)], "target": ["G"] * 5, "score": range(5)})
    api = excel_api(tmp_path, monkeypatch, table, report_name="rna_report")
    result = api.post_process_probes(table, raw_csv=True)
    files = api.output_files
    assert files["raw_xlsx"] is None and files["xlsx"] is None and not list(tmp_path.glob("*.xlsx"))
    for name in (files["raw_tsv"], files["tsv"]):
        assert pd.read_csv(tmp_path / name, sep="\t")["probe_id"].tolist() == [f"p{i}" for i in range(5)]
    assert [(e["table"], e["rows"], e["xlsx"]) for e in api.excel_limits] == [("raw", 5, None), ("final", 5, None)]
    assert len(result) == 5 and api.generate_report(result)["html_reports"]
