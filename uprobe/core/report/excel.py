"""Excel export and shared sequence/attribute color groups."""
from io import BytesIO
import json
import math
from pathlib import Path

import pandas as pd
import xlsxwriter

from .formatting import prepare_result_table, sequence_groups, column_group


PALETTE = [
    ("E5EDFA", "F3F6FC", "C9D8F0"),
    ("E1F0EA", "F1F8F4", "BFDCCC"),
    ("FCEEDB", "FEF8EE", "EBD4AF"),
    ("EEE5F7", "F7F3FC", "D9C8EA"),
    ("F8E4E9", "FCF3F5", "EAC8D2"),
    ("DDF1F4", "F0F9FA", "B9DDE2"),
]


def column_styles(columns, protocol):
    groups = sequence_groups(protocol)
    styles = {}
    for column in columns:
        group = column_group(column, protocol)
        if group is None:
            continue
        sequence, attribute, header = PALETTE[groups.index(group) % len(PALETTE)]
        styles[column] = {"group": group, "background": sequence if column == group else attribute,
                          "header": header, "sequence": column == group or column.startswith(group + ".")}
    return styles


def workbook_bytes(df, protocol):
    table = prepare_result_table(df, protocol)
    stream = BytesIO()
    book = xlsxwriter.Workbook(stream, {"in_memory": True, "strings_to_formulas": False,
                                       "strings_to_urls": False})
    sheet = book.add_worksheet("Probes")
    sheet.freeze_panes(1, min(2, len(table.columns)))
    sheet.set_zoom(90)
    sheet.hide_gridlines(2)
    styles = column_styles(table.columns, protocol)
    formats = {}
    for index, column in enumerate(table.columns):
        style = styles.get(column, {})
        body = {"font_name": "Consolas" if style.get("sequence") else "Calibri",
                "font_size": 11, "valign": "vcenter", "border": 1, "border_color": "#BBC7D2"}
        if style:
            body["bg_color"] = "#" + style["background"]
        formats[column] = book.add_format(body)
        header = book.add_format({"bold": True, "font_color": "#25313C", "font_size": 11,
                                  "bg_color": "#" + style.get("header", "E9EEF3"),
                                  "border": 1, "border_color": "#BBC7D2"})
        sheet.write_string(0, index, str(column), header)
        sheet.set_column(index, index, 48 if style.get("sequence") else 27 if column in {"target", "probe_id", "transcript_names"} else 20)
    sheet.set_row(0, 26)
    for row_index, row in enumerate(table.itertuples(index=False, name=None), 1):
        for column_index, (column, value) in enumerate(zip(table.columns, row)):
            if isinstance(value, (list, tuple, dict, set)):
                value = json.dumps(list(value) if isinstance(value, set) else value, default=str)
            if value is None or pd.isna(value) or isinstance(value, float) and not math.isfinite(value):
                sheet.write_blank(row_index, column_index, None, formats[column])
            elif isinstance(value, str):
                if len(value) > 32767:
                    raise ValueError(f"Excel cell exceeds 32767 characters: {column}")
                sheet.write_string(row_index, column_index, value, formats[column])
            else:
                sheet.write(row_index, column_index, value, formats[column])
    if len(table.columns):
        sheet.autofilter(0, 0, len(table), len(table.columns) - 1)
    book.close()
    return stream.getvalue()


def save_xlsx(df, protocol, path):
    path = Path(path)
    path.write_bytes(workbook_bytes(df, protocol))
    return path
