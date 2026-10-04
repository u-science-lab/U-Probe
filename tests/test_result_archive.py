import importlib.util
from pathlib import Path
import tempfile
import unittest
import zipfile

spec = importlib.util.spec_from_file_location(
    "result_archive", Path(__file__).parents[1] / "uprobe/http/utils/result_archive.py"
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class ResultArchiveTests(unittest.TestCase):
    def test_legacy_archive_keeps_report_and_both_tables_intact(self):
        with tempfile.TemporaryDirectory() as directory:
            archive = Path(directory) / "results.zip"
            results = {"report.html": b"<html>report</html>",
                       "probes_raw.csv": b"sequence\nACTG\n",
                       "probes.csv": b"sequence\nACTG\n"}
            with zipfile.ZipFile(archive, "w") as output:
                for name, content in results.items():
                    output.writestr(name, content)
                for name in ("protocol.yaml", "merged_genomes.yaml", "run.log", "data.fa"):
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
