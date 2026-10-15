import functools
import gzip
import http.server
import json
import random
import shutil
import subprocess
import threading
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml
from click.testing import CliRunner

from uprobe.core import genome_fetch as gf
from uprobe.core.genome_fetch import GenomeFetchError

SPECIES_HEADER = ("#name\tspecies\tdivision\ttaxonomy_id\tassembly\tassembly_accession\tgenebuild\t"
                  "variation\tmicroarray\tpan_compara\tpeptide_compara\tgenome_alignments\t"
                  "other_alignments\tcore_db\tspecies_id\n")


def _species_row(name, species, division, taxid, assembly, core_db):
    return f"{name}\t{species}\t{division}\t{taxid}\t{assembly}\tGCA_1\tbuild\tN\tN\tN\tN\tN\tN\t{core_db}\t1\t\n"


def _reference_sum(data: bytes):
    # BSD sum as in coreutils sum.c (bsd_sum_stream), independent of the module.
    checksum = 0
    for byte in data:
        checksum = (checksum >> 1) + ((checksum & 1) << 15)
        checksum = (checksum + byte) & 0xFFFF
    return checksum, (len(data) + 1023) // 1024


def _fasta(seed, n=20000):
    rng = random.Random(seed)
    seq = "".join(rng.choice("ACGT") for _ in range(n))
    return ">1 dna:chromosome\n" + "\n".join(seq[i:i + 60] for i in range(0, n, 60)) + "\n"


GTF = '1\tensembl\texon\t1\t100\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "GENE1";\n'


def _write_dir(directory: Path, files: dict):
    """Write gzipped files plus a CHECKSUMS file in Ensembl's format."""
    directory.mkdir(parents=True, exist_ok=True)
    lines = []
    for name, text in files.items():
        data = gzip.compress(text.encode(), mtime=0)
        (directory / name).write_bytes(data)
        checksum, blocks = _reference_sum(data)
        lines.append(f"{checksum:05d} {blocks} {name}\n")
    (directory / "README").write_text("readme\n")
    lines.append("{:05d} {} README\n".format(*_reference_sum(b"readme\n")))
    (directory / "CHECKSUMS").write_text("".join(lines))


class _RangeHandler(http.server.SimpleHTTPRequestHandler):
    """Static files with single 'bytes=N-' Range support, recording requests."""

    def log_message(self, *args):
        pass

    def do_GET(self):
        self.server.requests.append((self.path, self.headers.get("Range")))
        path = Path(self.translate_path(self.path))
        if not path.is_file():
            self.send_error(404)
            return
        data, start = path.read_bytes(), 0
        if self.headers.get("Range"):
            start = int(self.headers["Range"].split("=")[1].split("-")[0])
            if start >= len(data):
                self.send_error(416)
                return
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {start}-{len(data) - 1}/{len(data)}")
        else:
            self.send_response(200)
        self.send_header("Content-Length", str(len(data) - start))
        self.end_headers()
        self.wfile.write(data[start:])


@pytest.fixture
def ensembl(tmp_path):
    root = tmp_path / "server"
    ens = root / "ens" / "pub"
    (ens / "release-116").mkdir(parents=True)
    (ens / "VERSION").write_text("116\n")
    (ens / "release-116" / "species_EnsemblVertebrates.txt").write_text(
        SPECIES_HEADER
        + _species_row("Human", "homo_sapiens", "EnsemblVertebrates", 9606, "GRCh38.p14", "homo_sapiens_core_116_38")
        + _species_row("Caenorhabditis elegans (Nematode, N2)", "caenorhabditis_elegans", "EnsemblVertebrates",
                       6239, "WBcel235", "caenorhabditis_elegans_core_116_282")
        + _species_row("Mouse", "mus_musculus", "EnsemblVertebrates", 10090, "GRCm39", "mus_musculus_core_116_39"))
    human = SimpleNamespace(primary=_fasta(1), toplevel=_fasta(2), gtf=GTF)
    _write_dir(ens / "release-116/fasta/homo_sapiens/dna", {
        "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz": human.primary,
        "Homo_sapiens.GRCh38.dna.toplevel.fa.gz": human.toplevel,
        "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz": human.primary.lower(),
        "Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz": human.primary,
    })
    _write_dir(ens / "release-116/gtf/homo_sapiens", {
        "Homo_sapiens.GRCh38.116.gtf.gz": human.gtf,
        "Homo_sapiens.GRCh38.116.chr.gtf.gz": "chr only\n",
        "Homo_sapiens.GRCh38.116.abinitio.gtf.gz": "abinitio\n",
    })
    worm = SimpleNamespace(toplevel=_fasta(3, 3000), gtf=GTF)
    _write_dir(ens / "release-116/fasta/caenorhabditis_elegans/dna", {
        "Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz": worm.toplevel,
        "Caenorhabditis_elegans.WBcel235.dna_rm.toplevel.fa.gz": worm.toplevel.lower(),
    })
    # Worm/fly/yeast GTFs in Ensembl carry the Ensembl Genomes release number.
    _write_dir(ens / "release-116/gtf/caenorhabditis_elegans", {
        "Caenorhabditis_elegans.WBcel235.63.gtf.gz": worm.gtf,
    })
    eg = root / "eg" / "pub"
    (eg / "fungi/release-63").mkdir(parents=True)
    (eg / "VERSION").write_text("63\n")
    (eg / "fungi/release-63/species_EnsemblFungi.txt").write_text(
        SPECIES_HEADER + _species_row("Fungus x", "fungus_x_gca_1", "EnsemblFungi", 42, "ASM1v1",
                                      "fungi_ascomycota1_collection_core_63_116_1"))
    collection = "fungi_ascomycota1_collection/fungus_x_gca_1"
    _write_dir(eg / f"fungi/release-63/fasta/{collection}/dna",
               {"Fungus_x_gca_1.ASM1v1.dna.toplevel.fa.gz": _fasta(4, 2000)})
    _write_dir(eg / f"fungi/release-63/gtf/{collection}", {"Fungus_x_gca_1.ASM1v1.63.gtf.gz": GTF})
    (root / "rest/info/genomes").mkdir(parents=True)
    (root / "rest/info/data").write_text(json.dumps({"releases": [117]}))
    (root / "rest/info/genomes/arabidopsis_thaliana").write_text(
        json.dumps({"name": "arabidopsis_thaliana", "division": "EnsemblPlants"}))

    server = http.server.ThreadingHTTPServer(
        ("127.0.0.1", 0), functools.partial(_RangeHandler, directory=str(root)))
    server.requests = []
    thread = threading.Thread(target=server.serve_forever, kwargs={"poll_interval": 0.05}, daemon=True)
    thread.start()
    url = f"http://127.0.0.1:{server.server_address[1]}"
    yield SimpleNamespace(root=root, ens=f"{url}/ens/pub", eg=f"{url}/eg/pub", rest=f"{url}/rest",
                          url=url, requests=server.requests, human=human, worm=worm)
    server.shutdown()
    server.server_close()


@pytest.fixture
def fetch(ensembl, tmp_path):
    def run(species="homo_sapiens", **kwargs):
        kwargs.setdefault("mirror", ensembl.ens)
        kwargs.setdefault("rest_url", None)
        kwargs.setdefault("jellyfish", False)
        return gf.fetch_genome(species, out_dir=tmp_path / "genomes", genomes_yaml=tmp_path / "genomes.yaml",
                               **kwargs)
    return run


def _downloads(ensembl):
    return [path for path, _ in ensembl.requests if path.endswith(".gz")]


def test_bsd_sum_matches_reference_values(tmp_path, monkeypatch):
    sample = tmp_path / "sample.bin"
    sample.write_bytes(bytes(range(256)) * 20 + b"U-Probe")
    empty = tmp_path / "empty.bin"
    empty.write_bytes(b"")
    monkeypatch.setattr(gf.shutil, "which", lambda name: None)
    # Expected values were produced by the system `sum` (BSD algorithm).
    assert gf.bsd_sum(sample) == (48407, 6)
    assert gf.bsd_sum(empty) == (0, 0)


@pytest.mark.skipif(shutil.which("sum") is None, reason="no system sum")
def test_bsd_sum_system_and_python_paths_agree(tmp_path):
    big = tmp_path / "big.bin"
    big.write_bytes(random.Random(0).randbytes(gf._PY_SUM_LIMIT + 4321))
    out = subprocess.run(["sum", str(big)], capture_output=True, text=True, check=True).stdout.split()
    assert gf.bsd_sum(big) == gf._bsd_sum_python(big) == (int(out[0]), int(out[1]))


def test_parse_checksums_and_select_files():
    sums = gf.parse_checksums("22450 861294 A.GRCh38.dna.primary_assembly.fa.gz\n"
                              "01120 914531 A.GRCh38.dna.toplevel.fa.gz\n"
                              "18525 468971 A.GRCh38.dna_rm.primary_assembly.fa.gz\n\nbroken line\n")
    assert sums["A.GRCh38.dna.primary_assembly.fa.gz"] == (22450, 861294)
    assert len(sums) == 3
    gtfs = ["A.GRCh38.116.chr.gtf.gz", "A.GRCh38.116.gtf.gz", "A.GRCh38.116.abinitio.gtf.gz", "README"]
    assert gf.select_files(sums, gtfs, 116) == (
        "primary_assembly", "A.GRCh38.dna.primary_assembly.fa.gz", "A.GRCh38.116.gtf.gz", "GRCh38")
    toplevel_only = ["B.BDGP6.54.dna.toplevel.fa.gz", "B.BDGP6.54.dna_sm.toplevel.fa.gz"]
    assert gf.select_files(toplevel_only, ["B.BDGP6.54.63.chr.gtf.gz", "B.BDGP6.54.63.gtf.gz"], 116) == (
        "toplevel", "B.BDGP6.54.dna.toplevel.fa.gz", "B.BDGP6.54.63.gtf.gz", "BDGP6.54")
    with pytest.raises(GenomeFetchError, match="toplevel"):
        gf.select_files(["B.x.dna_sm.toplevel.fa.gz"], [], 116)


def test_release_resolution(ensembl):
    assert gf.resolve_release("latest", mirror=ensembl.ens, rest_url=None) == 116
    assert gf.resolve_release("latest", "fungi", mirror=ensembl.eg, rest_url=None) == 63
    assert gf.resolve_release("115") == 115
    assert gf.resolve_release("release-110") == 110
    with pytest.raises(GenomeFetchError, match="Invalid release"):
        gf.resolve_release("newest")
    # Without a VERSION file, the REST /info/data answer is used.
    assert gf.resolve_release("latest", mirror=f"{ensembl.url}/missing", rest_url=ensembl.rest) == 117
    with pytest.raises(GenomeFetchError, match="explicit --release"):
        gf.resolve_release("latest", mirror=f"{ensembl.url}/missing", rest_url=None)


def test_species_resolution_and_suggestions(ensembl):
    for query in ("homo_sapiens", "Human", "Homo sapiens", "9606"):
        plan = gf.plan_fetch(query, mirror=ensembl.ens, rest_url=None)
        assert (plan.species, plan.release, plan.assembly) == ("homo_sapiens", 116, "GRCh38")
    with pytest.raises(GenomeFetchError, match=r"Did you mean: homo_sapiens \(Human\)"):
        gf.plan_fetch("homo_sapien", mirror=ensembl.ens, rest_url=None)
    with pytest.raises(GenomeFetchError, match="--division plants --species arabidopsis_thaliana"):
        gf.plan_fetch("arabidopsis_thaliana", mirror=ensembl.ens, rest_url=ensembl.rest)
    with pytest.raises(GenomeFetchError, match="CHECKSUMS not found"):
        gf.plan_fetch("mus_musculus", mirror=ensembl.ens, rest_url=None)
    rows = gf.search_species("mouse", gf.load_species_table("vertebrates", 116, ensembl.ens))
    assert rows[0]["species"] == "mus_musculus"


def test_fetch_registers_primary_assembly(ensembl, fetch, tmp_path):
    result = fetch(release="latest")
    name = result["name"]
    assert name == "homo_sapiens_GRCh38_e116"
    directory = tmp_path / "genomes" / name
    entry = yaml.safe_load((tmp_path / "genomes.yaml").read_text())[name]
    assert entry["fasta"] == str(directory.resolve() / "Homo_sapiens.GRCh38.dna.primary_assembly.fa")
    assert Path(entry["fasta"]).read_text() == ensembl.human.primary
    assert Path(entry["gtf"]).read_text() == ensembl.human.gtf
    assert entry["align_index"] == ["bowtie2"] and entry["jellyfish"] is False
    assert entry["source"]["release"] == 116
    assert entry["source"]["sequence"] == "primary_assembly"
    assert entry["source"]["fasta_url"].endswith(
        "/release-116/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")
    assert sorted(p.name for p in directory.iterdir()) == [
        "Homo_sapiens.GRCh38.116.gtf", "Homo_sapiens.GRCh38.dna.primary_assembly.fa", gf.MANIFEST_NAME]
    manifest = json.loads((directory / gf.MANIFEST_NAME).read_text())
    assert manifest["release"] == 116
    assert manifest["files"]["gtf"]["file"] == "Homo_sapiens.GRCh38.116.gtf.gz"
    # The entry is in the format UProbeAPI expects.
    from uprobe.core.api import UProbeAPI
    api = UProbeAPI({"genome": name}, tmp_path / "genomes.yaml", tmp_path / "results")
    assert api.genome["fasta"] == entry["fasta"]


def test_toplevel_fallback(ensembl, fetch):
    result = fetch("caenorhabditis_elegans", release=116)
    entry = result["entry"]
    assert result["plan"].sequence == "toplevel"
    assert entry["fasta"].endswith("Caenorhabditis_elegans.WBcel235.dna.toplevel.fa")
    assert entry["gtf"].endswith("Caenorhabditis_elegans.WBcel235.63.gtf")
    assert Path(entry["fasta"]).read_text() == ensembl.worm.toplevel


def test_ensembl_genomes_collection_species(ensembl, fetch):
    result = fetch("fungus_x_gca_1", division="fungi", mirror=ensembl.eg)
    assert result["name"] == "fungus_x_gca_1_ASM1v1_eg63"
    assert "/fungi/release-63/fasta/fungi_ascomycota1_collection/fungus_x_gca_1/dna/" in result["entry"]["source"]["fasta_url"]
    assert Path(result["entry"]["gtf"]).read_text() == GTF


def test_corrupted_download_is_rejected(ensembl, fetch, tmp_path):
    served = ensembl.root / "ens/pub/release-116/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    served.write_bytes(gzip.compress(b">1\nACGTTTTT\n", mtime=0))
    with pytest.raises(GenomeFetchError, match="Checksum mismatch"):
        fetch()
    directory = tmp_path / "genomes" / "homo_sapiens_GRCh38_e116"
    assert not any(p.name.startswith("Homo_sapiens.GRCh38.dna") for p in directory.iterdir())
    assert not (tmp_path / "genomes.yaml").exists()
    assert len([p for p in _downloads(ensembl) if "primary_assembly" in p]) == 2  # retried once


def test_rerun_skips_verified_files(ensembl, fetch, tmp_path):
    fetch()
    text = (tmp_path / "genomes.yaml").read_text()
    ensembl.requests.clear()
    result = fetch()
    assert _downloads(ensembl) == []
    assert result["changed"] is False
    assert (tmp_path / "genomes.yaml").read_text() == text
    # A damaged local file (size changed) is fetched again.
    Path(result["entry"]["gtf"]).write_text("truncated")
    fetch()
    assert [p.rsplit("/", 1)[1] for p in _downloads(ensembl)] == ["Homo_sapiens.GRCh38.116.gtf.gz"]
    assert Path(result["entry"]["gtf"]).read_text() == GTF


def test_partial_download_is_resumed(ensembl, fetch, tmp_path):
    name = "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    data = (ensembl.root / "ens/pub/release-116/fasta/homo_sapiens/dna" / name).read_bytes()
    directory = tmp_path / "genomes" / "homo_sapiens_GRCh38_e116"
    directory.mkdir(parents=True)
    (directory / f"{name}.part").write_bytes(data[:1000])
    result = fetch()
    assert (f"/ens/pub/release-116/fasta/homo_sapiens/dna/{name}", "bytes=1000-") in ensembl.requests
    assert Path(result["entry"]["fasta"]).read_text() == ensembl.human.primary
    assert not list(directory.glob("*.part"))


def test_genomes_yaml_merge_keeps_other_entries_and_settings(fetch, tmp_path):
    path = tmp_path / "genomes.yaml"
    path.write_text("# my genomes\ncustom:\n  fasta: /data/custom.fa\n  gtf: /data/custom.gtf\n"
                    "  align_index:\n  - bowtie2\n  jellyfish: false\n")
    name = fetch()["name"]
    text = path.read_text()
    assert text.startswith("# my genomes\ncustom:")
    data = yaml.safe_load(text)
    assert list(data) == ["custom", name]
    # User settings of an existing entry survive an update.
    data[name].update({"align_index": ["bowtie2", "blast"], "jellyfish": True, "notes": "mine"})
    data[name]["description"] = "old description"
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    result = fetch(jellyfish=None)
    entry = yaml.safe_load(path.read_text())[name]
    assert result["changed"] is True
    assert entry["align_index"] == ["bowtie2", "blast"]
    assert entry["jellyfish"] is True and entry["notes"] == "mine"
    assert entry["description"] == "Human GRCh38 (Ensembl release 116, primary_assembly)"
    assert yaml.safe_load(path.read_text())["custom"]["fasta"] == "/data/custom.fa"


def test_directory_holding_other_release_needs_force(fetch, tmp_path):
    result = fetch(name="hg")
    directory = result["directory"]
    manifest = json.loads((directory / gf.MANIFEST_NAME).read_text())
    manifest["release"] = 115
    (directory / gf.MANIFEST_NAME).write_text(json.dumps(manifest))
    (directory / "bowtie2_genome").mkdir()
    (directory / "bowtie2_genome" / "x.1.bt2").write_text("stale")
    with pytest.raises(GenomeFetchError, match="already holds homo_sapiens GRCh38 .*release 115"):
        fetch(name="hg")
    fetch(name="hg", dry_run=True, force=True)
    assert (directory / "bowtie2_genome").exists()  # a dry run deletes nothing
    fetch(name="hg", force=True)
    assert not (directory / "bowtie2_genome").exists()
    assert json.loads((directory / gf.MANIFEST_NAME).read_text())["release"] == 116


def test_build_index_reuses_build_genome(fetch, monkeypatch):
    from unittest.mock import Mock
    import uprobe.core.api as api
    build = Mock()
    monkeypatch.setattr(api, "build_genome", build)
    result = fetch(build_index=True, threads=3)
    build.assert_called_once()
    assert build.call_args.args[0]["fasta"] == result["entry"]["fasta"]
    assert build.call_args.kwargs == {"threads": 3}
    build.side_effect = RuntimeError("bowtie2-build: not found")
    with pytest.raises(GenomeFetchError, match="downloaded and registered, but building indices failed"):
        fetch(build_index=True)


def test_cli_fetch_list_and_search(ensembl, tmp_path):
    from uprobe.core.cli import cli
    runner = CliRunner()
    genomes = tmp_path / "cli_genomes.yaml"
    common = ["--species", "human", "--mirror", ensembl.ens, "--out", str(tmp_path / "g"), "--genomes", str(genomes)]
    dry = runner.invoke(cli, ["genome", "fetch", *common, "--dry-run"])
    assert dry.exit_code == 0, dry.output
    assert yaml.safe_load(dry.output)["plan"]["release"] == 116
    assert not genomes.exists() and _downloads(ensembl) == []

    done = runner.invoke(cli, ["genome", "fetch", *common, "--no-jellyfish", "--align-index", "bowtie2,blast"])
    assert done.exit_code == 0, done.output
    entry = yaml.safe_load(done.output)["homo_sapiens_GRCh38_e116"]
    assert entry["align_index"] == ["bowtie2", "blast"]

    listed = runner.invoke(cli, ["genome", "list", "-g", str(genomes), "--json"])
    assert listed.exit_code == 0, listed.output
    [row] = json.loads(listed.output)
    assert row["fasta_exists"] and row["gtf_exists"]
    assert row["indexes"]["bowtie2"] == {"genome": False, "transcript": False}
    assert row["source"]["release"] == 116
    text = runner.invoke(cli, ["genome", "list", "-g", str(genomes)])
    assert "bowtie2 index: genome missing, transcript missing" in text.output

    found = runner.invoke(cli, ["genome", "search", "9606", "--mirror", ensembl.ens])
    assert found.exit_code == 0 and found.output.startswith("homo_sapiens\tHuman")
    bad = runner.invoke(cli, ["genome", "fetch", *common, "--release", "abc"])
    assert bad.exit_code == 1
