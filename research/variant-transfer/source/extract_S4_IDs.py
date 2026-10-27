#!/usr/bin/env python3
"""Extract only S4's text Variant IDs, never its numeric or categorical outcomes."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as E
import zipfile

S = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
EXPECTED = "723c4a0e1cfdb2c9c6ae55b3fa9206855dd95e0241319b9fc8d0831d3b42f9c2"


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def run(source, out):
    if digest(source) != EXPECTED:
        raise ValueError("S4 source hash mismatch")
    rows = []
    wanted = set()
    blank_rows = []
    with zipfile.ZipFile(source) as z:
        with z.open("xl/worksheets/sheet1.xml") as h:
            for _, row in E.iterparse(h, events=("end",)):
                if row.tag != S + "row":
                    continue
                nr = int(row.attrib["r"])
                if nr <= 2:
                    row.clear()
                    continue
                cells = [c for c in row.findall(S + "c") if c.attrib["r"] == f"A{nr}"]
                if len(cells) != 1:
                    if any(c.find(S + "v") is not None or c.find(S + "is") is not None
                           for c in row.findall(S + "c")):
                        raise ValueError(f"populated row lacks one authoritative Variant ID at row{nr}")
                    blank_rows.append(nr)
                    row.clear()
                    continue
                cell = cells[0]
                if cell.attrib.get("t") != "s" or cell.find(S + "f") is not None:
                    raise ValueError(f"Variant ID is not ordinary shared-string metadata at row{nr}")
                index = int(cell.findtext(S + "v"))
                wanted.add(index)
                rows.append({"source_row_1based": nr, "shared_string_index": index})
                row.clear()
        strings = {}
        with z.open("xl/sharedStrings.xml") as h:
            i = 0
            for _, node in E.iterparse(h, events=("end",)):
                if node.tag != S + "si":
                    continue
                if i in wanted:
                    strings[i] = "".join(n.text or "" for n in node.findall(".//" + S + "t"))
                i += 1
                node.clear()
                if wanted <= strings.keys():
                    break
    if wanted != strings.keys():
        raise ValueError("selected Variant ID shared string missing")
    for row in rows:
        row["id"] = strings[row.pop("shared_string_index")]
        if not row["id"] or row["id"] != row["id"].strip():
            raise ValueError("empty ID or ID requiring silent whitespace repair")
    counts = Counter(r["id"] for r in rows)
    result = {
        "schema_version": 1, "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "source_sha256": EXPECTED, "scope": "Column A text Variant metadata only; no other cell value read",
        "sheet": "S4.MPRA results", "source_member": "xl/worksheets/sheet1.xml",
        "data_row_count": len(rows), "unique_ID_count": len(counts),
        "duplicate_IDs": {k: v for k, v in counts.items() if v > 1},
        "blank_format_rows": blank_rows, "rows": rows,
        "ordered_ID_sha256": hashlib.sha256("\n".join(r["id"] for r in rows).encode()).hexdigest(),
        "implementation_sha256": digest(__file__), "numeric_or_outcome_column_values_read": 0,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.GzipFile(filename=str(out), mode="xb", mtime=0) as f:
        f.write((json.dumps(result, separators=(",", ":")) + "\n").encode())
    print(json.dumps({k: v for k, v in result.items() if k != "rows"}, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    run(args.source, args.out)


if __name__ == "__main__":
    main()
