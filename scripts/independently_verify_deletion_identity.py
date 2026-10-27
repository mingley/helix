#!/usr/bin/env python3
"""Separate Decimal/sequence implementation of the HBB deletion-identity audit."""
import collections
import csv
from decimal import Decimal
import hashlib
import json
from pathlib import Path


def compute(root):
    table = root / "research/forensic-sources/GRCh38_HBB.tsv"
    report = json.loads((root / "research/regulatory-forensics.json").read_text())
    raw = table.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == report["input_commitments"]["original_osf_table_sha256"]
    sequence = "".join((root / "datasets/hbb-promoter-assay/target.fasta").read_text().splitlines()[1:])
    assert len(sequence) == 187
    all_groups = collections.defaultdict(list)
    qc_groups = collections.defaultdict(list)
    for row in csv.DictReader(raw.decode().splitlines(), delimiter="\t"):
        if row["Alt"] != "-":
            continue
        offset = int(row["Position"]) - 1 - 5227021
        assert sequence[offset] == row["Ref"]
        altered = sequence[:offset] + sequence[offset + 1:]
        assert len(altered) == 186
        item = dict(position=offset + 1, score=Decimal(row["Value"]), tags=int(row["Tags"]))
        all_groups[altered].append(item)
        if item["tags"] >= 10 and Decimal(row["DNA"]) > 0 and Decimal(row["RNA"]) > 0:
            qc_groups[altered].append(item)

    def summary(groups):
        n = sum(map(len, groups.values()))
        total = Decimal(0)
        duplicate_pairs = []
        for rows in groups.values():
            center = sum((r["score"] for r in rows), Decimal(0)) / len(rows)
            total += sum(((r["score"] - center) ** 2 for r in rows), Decimal(0))
            if len(rows) > 1:
                duplicate_pairs.append(dict(
                    positions=[r["position"] for r in rows], tags=[r["tags"] for r in rows],
                    one_below_ten_one_at_least_ten=sorted(r["tags"] >= 10 for r in rows) == [False, True]))
        return dict(rows=n, distinct_altered_sequences=len(groups), sum_squared_deviations=str(total),
                    minimum_observed_label_mse=str(total / n), duplicate_groups=duplicate_pairs)

    raw_result, qc_result = summary(all_groups), summary(qc_groups)
    assert raw_result["rows"] == 62 and raw_result["distinct_altered_sequences"] == 53
    assert Decimal(raw_result["sum_squared_deviations"]) == Decimal("4.480450")
    assert all(group["one_below_ten_one_at_least_ten"] for group in raw_result["duplicate_groups"])
    assert qc_result["rows"] == qc_result["distinct_altered_sequences"] == 35
    assert Decimal(qc_result["minimum_observed_label_mse"]) == 0
    return dict(schema_version=1, method="Independent original-TSV/whole-edited-sequence grouping with exact Decimal arithmetic; no shared reconciliation code",
                raw=raw_result, predefined_quality_pass=qc_result,
                scope="Finite source-estimate identity diagnostic; no biological irreducibility or independent replication claim")


if __name__ == "__main__":
    print(json.dumps(compute(Path(__file__).resolve().parents[1]), indent=2))
