#!/usr/bin/env python3
"""Normalize published MERFISH workbooks; never generate target-binding candidates."""

import argparse
import csv
import hashlib
import itertools
import json
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import yaml

SOURCE_SHA256 = {
    "elife-90029-supp1.xlsx": "bdabf93e017f146b8ad50c4dd724bf23d5f093806e894a0ac64814bd2f6d28b6",
    "elife-90029-supp2.xlsx": "e259948fca6b4184f51656936898a9c994458b8e1ebb94ce127ab85c67a12d9a",
    "elife-90029-supp3.xlsx": "105f725d58ddb1b4c81d288fe014ee2cc35fc04ba6468a5473b3a7507929d612",
}
PREFIX = "CGCGGTGGTGGTACGACATGA"
SUFFIX = "AGCGCAACTAGGACGGACAA"
ROLE = "encoding_probe"


def read_workbook(path):
    """Read the single worksheet, supporting both strict and transitional OOXML."""
    with zipfile.ZipFile(path) as archive:
        namespace = ET.fromstring(archive.read("xl/workbook.xml")).tag.split("}")[0] + "}"
        strings = [
            "".join(node.itertext())
            for node in ET.fromstring(archive.read("xl/sharedStrings.xml"))
        ]
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        rows = []
        for row in sheet.findall(namespace + "sheetData/" + namespace + "row"):
            values = {}
            for cell in row:
                value_node = cell.find(namespace + "v")
                value = value_node.text if value_node is not None else ""
                column = cell.get("r").rstrip("0123456789")
                values[column] = strings[int(value)] if cell.get("t") == "s" else value
            rows.append(values)
        return rows


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def variant_expression(key):
    return {
        "expr": (
            "encoding[target]['variants']"
            "[probe_index % encoding[target]['variant_count']]"
            "['" + key + "']"
        )
    }


def build_design(targets, encoding):
    return {
        "name": "MERFISH_2021_242gene_unfiltered_candidate_recovery",
        "genome": "mm10",
        "targets": [target["target_id"] for target in targets],
        "encoding": encoding,
        "extracts": {
            "target_region": {
                "source": "mRNA", "length": 30, "step": 1, "transcripts": "all"
            }
        },
        "probes": {
            ROLE: {
                "template": PREFIX + "{left}{binding}A{right}" + SUFFIX,
                "parts": {
                    "left": variant_expression("left"),
                    "binding": {"expr": "rc(target_region)"},
                    "right": variant_expression("right"),
                },
            }
        },
        "attributes": {},
        "post_process": {},
        "recovery_execution": {
            "encoding_variants": "cartesian",
            "transcript_policy": (
                "Each annotated isoform separately through UProbeAPI; "
                "original transcript mapping audited separately."
            ),
            "candidate_policy": (
                "All ACGT 30-mer transcript windows crossed with all 12 ordered distinct "
                "readout pairs from the four on-bits. No reference binding or complete "
                "encoding sequences are generation inputs."
            ),
            "final_selection": "Not evaluated",
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for filename, expected in SOURCE_SHA256.items():
        observed = hashlib.sha256((args.raw_dir / filename).read_bytes()).hexdigest()
        if observed != expected:
            raise ValueError("Source SHA256 mismatch: " + filename)
    code_rows, probe_rows, readout_rows = [
        read_workbook(args.raw_dir / ("elife-90029-supp%d.xlsx" % index))
        for index in (1, 2, 3)
    ]
    # Both workbooks contain an Excel-converted March1 gene cell.
    # ENSMUST00000072482 maps unambiguously to March1 in mm10 GTF.
    for row in code_rows[1:]:
        if row["B"] == "ENSMUST00000072482":
            row["A"] = "March1"

    bit_columns = list(code_rows[0])[2:]
    codebook = [
        {
            "target_id": row["A"].upper(),
            "transcript_id": row["B"],
            "barcode": "".join(row[column] for column in bit_columns),
            **{code_rows[0][column]: row[column] for column in bit_columns},
        }
        for row in code_rows[1:]
    ]
    targets = [
        {
            "target_id": row["target_id"],
            "annotation_gene_name": next(
                source["A"] for source in code_rows[1:]
                if source["B"] == row["transcript_id"]
            ),
            "transcript_id": row["transcript_id"],
            "barcode": row["barcode"],
        }
        for row in codebook if not row["target_id"].startswith("BLANK-")
    ]
    assert len(targets) == 242
    targets_by_transcript = {row["transcript_id"]: row for row in targets}
    readouts = [
        {
            "bit": row["A"], "readout_id": row["B"],
            "readout_sequence": row["C"],
            "encoding_sequence": reverse_complement(row["C"]), "dye": row["D"],
        }
        for row in readout_rows[1:]
    ]
    readout_by_sequence = {row["encoding_sequence"]: row["readout_id"] for row in readouts}
    sequence_by_readout = {row["readout_id"]: row["encoding_sequence"] for row in readouts}
    encoding = {}
    for target in targets:
        row = next(row for row in codebook if row["target_id"] == target["target_id"])
        bits = [
            code_rows[0][column] for column in bit_columns
            if row[code_rows[0][column]] == "1"
        ]
        assert len(bits) == 4
        variants = [
            {"left": sequence_by_readout[left], "right": sequence_by_readout[right]}
            for left, right in itertools.permutations(bits, 2)
        ]
        encoding[target["target_id"]] = {
            "variants": variants, "variant_count": len(variants)
        }

    reference, corrections, pairs = [], [], {}
    for source_row, row in enumerate(probe_rows[1:], 2):
        target = targets_by_transcript[row["B"]]
        gene, sequence = target["target_id"], row["C"].upper()
        if row["A"].upper() != gene:
            corrections.append([source_row, row["A"], gene, row["B"]])
        assert len(sequence) == 112 and set(sequence) <= set("ACGT")
        assert sequence[:21] == PREFIX
        assert sequence[71] == "A" and sequence[92:] == SUFFIX
        left = readout_by_sequence[sequence[21:41]]
        right = readout_by_sequence[sequence[72:92]]
        assert any(
            variant["left"] == sequence[21:41] and variant["right"] == sequence[72:92]
            for variant in encoding[gene]["variants"]
        )
        pairs.setdefault(gene, set()).add((left, right))
        reference.append({
            "gene": gene, "pair_id": "", "oligo_type": ROLE,
            "oligo_id": "%s_encoding_%05d" % (gene, source_row - 1),
            "barcode": target["barcode"], "full_sequence": sequence,
            "dna_sequence": sequence, "binding_sequence": sequence[41:71],
            "binding_length_nt": 30, "transcript_id": row["B"],
            "readout_left": left, "readout_right": right,
            "source_table": "eLife 90029 Supplementary file 2", "source_row": source_row,
        })
    assert len(reference) == 22005 and len(pairs) == 242 and len(readouts) == 22

    for name, rows in [
        ("targets", targets), ("reference_oligos", reference),
        ("codebook", codebook), ("readouts", readouts),
    ]:
        write_csv(args.output_dir / (name + ".csv"), rows)
    design = build_design(targets, encoding)
    (args.output_dir / "design.yaml").write_text(
        yaml.safe_dump(design, sort_keys=False, width=100)
    )
    print(json.dumps({
        "targets": len(targets), "reference_oligos": len(reference),
        "readouts": len(readouts), "codebook_rows": len(codebook),
        "corrected_gene_cells": corrections,
        "pair_counts": sorted(set(map(len, pairs.values()))),
        "source_sha256": SOURCE_SHA256,
    }))


if __name__ == "__main__":
    main()
