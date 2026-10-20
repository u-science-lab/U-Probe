import warnings

import pytest
from pyfaidx import Fasta
from uprobe.core.gen.fun import read_gtf, validate_targets, get_exon_seq, extract_trans_seqs


@pytest.mark.parametrize("attributes, expected", [
    ('gene_name "Preferred"; gene "Fallback";', "Preferred"),
    ('gene "Fallback";', "Fallback"),
    ('gene_name ""; gene "Fallback";', "Fallback"),
    ('gene_name "   "; gene "Fallback";', "Fallback"),
    ('gene\t"Fallback";', "Fallback"),
])
def test_attribute_fallback(tmp_path, attributes, expected):
    # A literal tab is not valid inside the ninth column; test whitespace with spaces.
    attributes = attributes.replace("\t", "  ")
    gtf = tmp_path / "test.gtf"
    gtf.write_text('1\ttest\texon\t1\t80\t.\t+\t.\tgene_id "ID1"; transcript_id "TX1"; ' + attributes + "\n")
    original = gtf.read_bytes()
    assert read_gtf(gtf).gene_name.tolist() == [expected]
    assert validate_targets([expected, "Absent"], gtf, DTF_NAME_FIX=True) == (False, [expected], ["Absent"])
    assert gtf.read_bytes() == original


def test_mixed_annotation_and_sequence_extraction(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">1\n" + "ACGT" * 50 + "\n")
    gtf = tmp_path / "test.gtf"
    rows = [
        '1\ttest\texon\t1\t80\t.\t+\t.\tgene_id "ID1"; transcript_id "TX1"; gene "Legacy";',
        '1\ttest\texon\t101\t180\t.\t-\t.\tgene_id "ID2"; transcript_id "TX2"; gene_name "Modern"; gene "Ignored";',
        '1\ttest\tgene\t181\t200\t.\t+\t.\tgene_id "ID3"; other_gene "FalseMatch";',
    ]
    gtf.write_text("\n".join(rows) + "\n")
    original = gtf.read_bytes()
    parsed = read_gtf(gtf)
    assert parsed.gene_name.iloc[:2].tolist() == ["Legacy", "Modern"]
    assert parsed.gene_name.iloc[2:].isna().all()
    assert validate_targets(["Legacy", "Modern"], gtf, DTF_NAME_FIX=True) == (True, ["Legacy", "Modern"], [])
    seqs = get_exon_seq(["Legacy", "Modern"], str(fa), str(gtf))
    assert seqs["Legacy"][0][2] == "ACGT" * 20
    assert seqs["Modern"][0][2] == "ACGT" * 20
    transcript = tmp_path / "transcript.fa"
    extract_trans_seqs(str(gtf), str(fa), str(transcript))
    with Fasta(str(transcript)) as transcripts:
        assert set(transcripts.keys()) == {"ID1_TX1", "ID2_TX2"}
    assert gtf.read_bytes() == original


def test_read_gtf_reads_chromosome_names_as_text(tmp_path):
    # A large GTF with numeric (1) and other (X) chromosome names: pandas typed the
    # column per chunk and every run and validate log said "DtypeWarning: Columns (0)
    # have mixed types" (human Ensembl GTF).
    line = '{}\ttest\texon\t1\t80\t.\t+\t.\tgene_id "{}"; gene_name "{}";\n'
    gtf = tmp_path / "large.gtf"
    gtf.write_text(line.format(1, "ID1", "G1") * 100000 + line.format("X", "ID2", "G2"))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        df = read_gtf(gtf)
    assert [str(w.message) for w in caught if "mixed types" in str(w.message)] == []
    assert df.chr.iloc[0] == "1" and df.chr.iloc[-1] == "X"
