import importlib.util
from pathlib import Path
import tempfile
import unittest
import zipfile
import io

spec = importlib.util.spec_from_file_location(
    "result_archive", Path(__file__).parents[1] / "uprobe/http/utils/result_archive.py"
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class ResultArchiveTests(unittest.TestCase):
    def test_empty_table_removed_but_raw_candidates_and_report_retained(self):
        def workbook(rows):
            stream = io.BytesIO()
            with zipfile.ZipFile(stream, 'w') as book:
                book.writestr('xl/worksheets/sheet1.xml',
                    '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>'
                    '<row r="1"><c t="inlineStr"><is><t>probe_id</t></is></c></row>' + rows + '</sheetData></worksheet>')
            return stream.getvalue()
        with tempfile.TemporaryDirectory() as directory:
            archive = Path(directory) / 'results.zip'
            raw = workbook('<row r="2"><c t="inlineStr"><is><t>p1</t></is></c></row>')
            with zipfile.ZipFile(archive, 'w') as output:
                output.writestr('empty.xlsx', workbook(''))
                output.writestr('styled_empty.xlsx', workbook('<row r="2"><c s="1" /></row>'))
                output.writestr('probes_raw.xlsx', raw)
                output.writestr('report.html', '<html>report</html>')
            module.restrict_result_archive(archive)
            with zipfile.ZipFile(archive) as output:
                self.assertEqual(set(output.namelist()), {'probes_raw.xlsx', 'report.html'})
                self.assertEqual(output.read('probes_raw.xlsx'), raw)
            original = archive.read_bytes()
            module.restrict_result_archive(archive)
            self.assertEqual(archive.read_bytes(), original)

    def test_legacy_archive_keeps_report_and_both_tables_intact(self):
        with tempfile.TemporaryDirectory() as directory:
            archive = Path(directory) / "results.zip"
            results = {"report.html": b"<html>report</html>",
                       "probes_raw.xlsx": b"sequence\nACTG\n",
                       "probes.xlsx": b"sequence\nACTG\n"}
            with zipfile.ZipFile(archive, "w") as output:
                for name, content in results.items():
                    output.writestr(name, content)
                for name in ("protocol.yaml", "merged_genomes.yaml", "run.log", "data.fa", "legacy.csv"):
                    output.writestr(name, "excluded")
            module.restrict_result_archive(archive)
            with zipfile.ZipFile(archive) as output:
                self.assertEqual(set(output.namelist()), set(results))
                self.assertIsNone(output.testzip())
                for name, content in results.items():
                    self.assertEqual(output.read(name), content)
            original = archive.read_bytes()
            module.restrict_result_archive(archive)
            self.assertEqual(archive.read_bytes(), original)
            self.assertEqual(list(Path(directory).iterdir()), [archive])

    def test_invalid_archive_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            archive = Path(directory) / "results.zip"
            archive.write_bytes(b"invalid zip")
            with self.assertRaises(zipfile.BadZipFile):
                module.restrict_result_archive(archive)
            self.assertEqual(archive.read_bytes(), b"invalid zip")


if __name__ == "__main__":
    unittest.main()
