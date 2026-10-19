"""The genomes.yaml entry the HTTP genome router writes for a genome a user adds."""
import os
import tempfile

import pytest
import yaml

# The HTTP modules resolve data_dir at import time; keep them off the configured server paths.
os.environ.setdefault("UPROBE_SERVER_ROOT", tempfile.mkdtemp(prefix="uprobe-http-test-"))

from uprobe.http.routers import genome as router  # noqa: E402


@pytest.mark.parametrize("installed", [True, False])
def test_new_genome_gets_bowtie_when_bowtie_build_is_installed(tmp_path, monkeypatch, installed):
    # As `uprobe genome fetch` does: Bowtie 1 for the exhaustive off-target search of short probes.
    which = router.shutil.which
    monkeypatch.setattr(router.shutil, "which", lambda name, *args, **kwargs: (
        ("/opt/bin/bowtie-build" if installed else None) if name == "bowtie-build" else which(name, *args, **kwargs)))
    monkeypatch.setattr(router, "get_user_genomes_yaml", lambda username: tmp_path / username / "genomes.yaml")
    monkeypatch.setattr(router, "get_user_genomes_dir", lambda username: tmp_path / username)
    router.update_genomes_yaml("toy", "alice", "add")
    entry = yaml.safe_load((tmp_path / "alice" / "genomes.yaml").read_text())["toy"]
    assert entry["align_index"] == (["bowtie2", "blast", "bowtie"] if installed else ["bowtie2", "blast"])
    assert entry["jellyfish"] is False and entry["fasta"] == f"{tmp_path / 'alice' / 'toy'}/toy.fa"
    # An existing entry is left as it is.
    monkeypatch.setattr(router.shutil, "which", lambda name, *args, **kwargs: None)
    router.update_genomes_yaml("toy", "alice", "add")
    assert yaml.safe_load((tmp_path / "alice" / "genomes.yaml").read_text())["toy"] == entry
