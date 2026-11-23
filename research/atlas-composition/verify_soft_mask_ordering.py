#!/usr/bin/env python3
"""Check the Rust common-weight result at every exact breakpoint and open interval."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

LABELS = ("cds", "exon_outside_cds", "transcript_span_outside_exon", "outside_transcript_span")
NUCLEAR = tuple(f"chr{i}" for i in range(1, 23)) + ("chrX", "chrY")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def gc(counts):
    called = sum(counts[b] for b in "acgt")
    return Fraction(counts["g"] + counts["c"], called) if called else None


def difference(first, second):
    return None if first is None or second is None else first - second


def interpolate(endpoints, w):
    unmasked, masked = endpoints
    if w == 0:
        return unmasked
    if w == 1:
        return masked
    if unmasked is None or masked is None:
        return None
    return (1 - w) * unmasked + w * masked


def exact(value):
    if value is None:
        return None
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def verify(report, study):
    require(report["operation"] == "atlas_composition_by_soft_mask", "expected a stratified atlas")
    nuclear = [c for c in report["chromosomes"] if c["name"] != "chrM"]
    require(tuple(c["name"] for c in nuclear) == NUCLEAR, "nuclear chromosome set or order differs")
    require(study["operation"] == "soft_mask_gc_ordering" and tuple(study["class_order"]) == LABELS,
            "unexpected ordering study or class order")
    require(study["atlas"] == report["atlas"] and study["reference"] == report["reference"],
            "study and count report have different source identities")
    modes = ("combined", "unmasked", "masked", "chromosome_common_mask_fraction")
    comparisons = {name: [] for name in modes}
    gaps = []
    crossings = []
    weights = {}
    breakpoints = {Fraction(0), Fraction(1)}
    for chrom in nuclear:
        classes = [chrom["classes"][key] for key in LABELS]
        totals = [sum(group[stratum][base] for group in classes for base in "acgt")
                  for stratum in ("unmasked", "masked")]
        weight = Fraction(totals[1], sum(totals)) if sum(totals) else None
        weights[chrom["name"]] = exact(weight)
        endpoints = [(gc(c["unmasked"]), gc(c["masked"])) for c in classes]
        pooled = [gc({base: c["unmasked"][base] + c["masked"][base] for base in "acgt"})
                  for c in classes]
        rates = (pooled, [r[0] for r in endpoints], [r[1] for r in endpoints],
                 [interpolate(r, weight) if weight is not None else None for r in endpoints])
        for mode, values in zip(modes, rates, strict=True):
            differences = [difference(a, b) for a, b in zip(values, values[1:])]
            comparisons[mode].append((chrom["name"], differences))
        gap_pairs = zip(comparisons["unmasked"][-1][1], comparisons["masked"][-1][1], strict=True)
        for i, (u, m) in enumerate(gap_pairs):
            gaps.append((u, m))
            if u is not None and m is not None and u * m < 0:
                crossing = u / (u - m)
                breakpoints.add(crossing)
                crossings.append({"contrast": {"chromosome": chrom["name"], "first_class": LABELS[i],
                                               "second_class": LABELS[i + 1]},
                                  "masked_fraction": exact(crossing), "strict_above": m > 0})
    require(study["chromosome_masked_called_fraction"] == weights, "chromosome mask fractions differ")
    require(study["common_mask_fraction_range"]["interior_crossings"] == crossings,
            "interior crossings differ")

    for mode, rows in comparisons.items():
        result = study["modes"][mode]
        defined = sum(all(v is not None for v in vs) for _, vs in rows)
        strict = sum(all(v is not None and v > 0 for v in vs) for _, vs in rows)
        require((defined, strict) == (result["chromosomes_with_defined_order"],
                                      result["chromosomes_in_strict_order"]), f"{mode}: order count differs")
        require(len(result["adjacent"]) == len(LABELS) - 1, f"{mode}: unexpected contrast count")
        for i, observed in enumerate(result["adjacent"]):
            column = [(name, values[i]) for name, values in rows]
            expected = {"first_class": LABELS[i], "second_class": LABELS[i + 1],
                        "greater": sum(v is not None and v > 0 for _, v in column),
                        "tied": [c for c, v in column if v == 0],
                        "reversed": [c for c, v in column if v is not None and v < 0],
                        "undefined": [c for c, v in column if v is None]}
            require(observed == expected, f"{mode}/{LABELS[i]}: adjacent comparisons differ")

    # Adding the proposed bounds tests spurious cuts even if they coincide
    # with neither an actual crossing nor an endpoint.
    interval = study["common_mask_fraction_range"]
    lower = interval["lower_bound"]
    upper = interval["upper_bound"]

    def ratio(bound):
        value = bound["masked_fraction"]
        return Fraction(int(value["numerator"]), int(value["denominator"]))

    lo, hi = ratio(lower), ratio(upper)
    require(0 <= lo <= 1 and 0 <= hi <= 1, "reported bounds outside [0, 1]")
    breakpoints.update((lo, hi))
    cuts = sorted(breakpoints)
    midpoints = [(left + right) / 2 for left, right in zip(cuts, cuts[1:])]
    for w in cuts + midpoints:
        values = [interpolate(pair, w) for pair in gaps]
        observed = all(v is not None and v > 0 for v in values)
        predicted = (not interval["empty"]
                     and (w > lo or (w == lo and lower["inclusive"]))
                     and (w < hi or (w == hi and upper["inclusive"])))
        require(observed == predicted, f"reported interval disagrees with raw counts at w={w}")
    return len(gaps), len(cuts), len(midpoints)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("counts", type=Path)
    parser.add_argument("ordering", type=Path)
    args = parser.parse_args()
    contrasts, cuts, interiors = verify(json.loads(args.counts.read_bytes()),
                                        json.loads(args.ordering.read_bytes()))
    print(f"Verified {contrasts} GC contrasts at {cuts} boundary points and {interiors} open intervals.")
