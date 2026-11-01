#!/usr/bin/env python3
"""Prepare fixed design inputs and evaluation-only probe-body references."""
import argparse
import csv
import re
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    def read(stem):
        # Official CSV contents may have been downloaded with a .txt suffix.
        path = args.raw_dir / (stem + ".csv")
        if not path.exists():
            path = args.raw_dir / (stem + ".txt")
        with path.open(newline="", encoding="utf-8-sig") as handle:
            return list(csv.DictReader(handle))

    def write(name, rows, columns):
        with (args.output_dir / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)

    # Target identities and non-targeting inputs come only from codebook/readouts.
    codebook = [row for row in read("channel 1 barcodekey-onoff")
                if not row["gene"].startswith("control-")]
    assert len(codebook) == 1194
    assert len({row["gene"].casefold() for row in codebook}) == 1194
    readouts = [dict(readout_id=row["Newname"], sequence=row["Sequence"],
                     round=row["Barcoded Hyb"], pseudocolor=row["Pseudocolors"])
                for row in read("channel 1-mRNA-seqFISH-readout-information")]
    assert len(readouts) == 64
    assert all(re.fullmatch("[ACGT]{15}", row["sequence"]) for row in readouts)
    lookup = {(int(row["round"]), int(row["pseudocolor"])): row
              for row in readouts}
    assert len(lookup) == 64
    write("targets.csv", [{"gene": row["gene"]} for row in codebook], ["gene"])
    write("codebook.csv", codebook, ["gene", "Hyb1", "Hyb2", "Hyb3", "Hyb4"])
    write("readouts.csv", readouts, ["readout_id", "sequence", "round", "pseudocolor"])

    # Published primary/full sequences are read only to prepare evaluation truth.
    probes = read("channel 1-mRNA-seqFISH-probes-without-primers")
    assert len(probes) == 37893
    by_gene = {row["gene"]: row for row in codebook}
    assert {row["name"] for row in probes} == set(by_gene)
    references = []
    for source_row, row in enumerate(probes, 2):
        gene, sequence, binding = row["name"], row["sequence"], row["primary_sequence"]
        assert re.fullmatch("[ACGT]{110}", sequence)
        assert re.fullmatch("[ACGT]{35}", binding) and sequence[30:65] == binding
        reads = [lookup[(j, int(by_gene[gene]["Hyb" + str(j)]))]["sequence"]
                 for j in range(1, 5)]
        arrangements = []
        for offset in range(4):
            r = reads[offset:] + reads[:offset]
            arrangements.append(r[0] + r[1] + binding + r[0] + r[2] + r[3])
        assert sequence in arrangements
        references.append(dict(gene=gene, oligo_type="primary_probe_body",
                               oligo_id=f"probe_body_{source_row - 1:05d}",
                               full_sequence=sequence, binding_sequence=binding,
                               source_row=source_row))
    assert len({(row["gene"].casefold(), row["full_sequence"])
                for row in references}) == 37893
    write("reference_oligos.csv", references,
          ["gene", "oligo_type", "oligo_id", "full_sequence", "binding_sequence", "source_row"])
    print("Prepared 1194 targets, 64 readouts and 37893 reference probe bodies.")


if __name__ == "__main__":
    main()
