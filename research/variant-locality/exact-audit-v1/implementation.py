#!/usr/bin/env python3
"""Independent exact, outcome-free audit of RC-3mer/GC SNV features.

No repository imports, floating-point arithmetic, experimental outcomes or fits.
Features: counts in 32 reverse-complement-canonical 3mer bins plus GC count.
For a 200nt sequence their denominators are 198 and 200 respectively.
"""
import argparse
import collections
import csv
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path

ALPHABET = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def rc(word):
    return word.translate(COMPLEMENT)[::-1]


WORDS = ["".join(chars) for chars in itertools.product(ALPHABET, repeat=3)]
CANONICAL = sorted({min(word, rc(word)) for word in WORDS})
INDEX = {word: number for number, word in enumerate(CANONICAL)}
BIN = {word: INDEX[min(word, rc(word))] for word in WORDS}


def counts(sequence):
    result = [0] * 33
    for start in range(len(sequence) - 2):
        result[BIN[sequence[start:start + 3]]] += 1
    result[-1] = sequence.count("G") + sequence.count("C")
    return tuple(result)


def subtract(after, before):
    return tuple(a - b for a, b in zip(after, before))


def local_delta(reference, alternative):
    # Independently accumulate lost and gained windows rather than subtracting
    # two precomputed full count vectors.
    result = [0] * 33
    mutant = reference[:2] + alternative + reference[3:]
    for start in (0, 1, 2):
        result[BIN[reference[start:start + 3]]] -= 1
        result[BIN[mutant[start:start + 3]]] += 1
    result[-1] = int(alternative in "GC") - int(reference[2] in "GC")
    return tuple(result)


def signature(vector):
    return ",".join(map(str, vector))


def digest(path):
    data = path.read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def write_csv(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def pair_checks(distance):
    # Entire sequence count subtraction, not a special interaction formula.
    # All words spanning the two edited sites and two flanks on either side.
    length = distance + 5
    table = {}
    for chars in itertools.product(ALPHABET, repeat=length):
        word = "".join(chars)
        table[word] = counts(word)
    tested = 0
    zero = 0
    first_nonzero = None
    second = 2 + distance
    for reference, base_counts in table.items():
        for alt1 in ALPHABET:
            if alt1 == reference[2]:
                continue
            first_word = reference[:2] + alt1 + reference[3:]
            first_counts = table[first_word]
            for alt2 in ALPHABET:
                if alt2 == reference[second]:
                    continue
                second_word = reference[:second] + alt2 + reference[second + 1:]
                double_word = first_word[:second] + alt2 + first_word[second + 1:]
                interaction = tuple(a - b - c + d for a, b, c, d in zip(
                    table[double_word], first_counts, table[second_word], base_counts))
                tested += 1
                if not any(interaction):
                    zero += 1
                elif first_nonzero is None:
                    first_nonzero = {"reference": reference, "mutant": double_word,
                        "positions_zero_based": [2, second],
                        "interaction_vector": list(interaction)}
                assert interaction[-1] == 0
    if distance >= 3:
        assert zero == tested
    return {"distance_nt": distance, "reference_contexts": len(table),
        "direct_full_count_comparisons": tested, "zero_interactions": zero,
        "nonzero_interactions": tested - zero, "first_nonzero": first_nonzero}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit("Refusing to overwrite existing results")
    assert len(CANONICAL) == 32
    classes = collections.defaultdict(list)
    rows = []
    rc_pairs = set()
    orientation_checks = 0
    for chars in itertools.product(ALPHABET, repeat=5):
        reference = "".join(chars)
        before = counts(reference)
        for alternative in ALPHABET:
            if alternative == reference[2]:
                continue
            mutant = reference[:2] + alternative + reference[3:]
            delta = local_delta(reference, alternative)
            assert delta == subtract(counts(mutant), before)
            assert sum(delta[:32]) == 0
            # Every interior changed nucleotide appears in exactly three windows.
            assert sum(delta[i] * (word.count("G") + word.count("C"))
                       for i, word in enumerate(CANONICAL)) == 3 * delta[-1]
            assert delta == local_delta(rc(reference), alternative.translate(COMPLEMENT))
            orientation_checks += 1
            pair_key = min(reference + ">" + mutant, rc(reference) + ">" + rc(mutant))
            rc_pairs.add(pair_key)
            entry = {"reference": reference, "alternative": alternative,
                "mutant": mutant, "rc_pair": pair_key, "delta": delta}
            classes[delta].append(entry)
            rows.append(entry)
    zero_vector = (0,) * 33
    zeros = classes[zero_vector]
    assert all(entry["mutant"] == rc(entry["reference"]) for entry in zeros)
    assert sum(entry["mutant"] == rc(entry["reference"]) for entry in rows) == len(zeros)
    outside_checks = 0
    # Exhaust every two-base flank on either side of each edited 5mer. All
    # full 9mer counts are computed separately from the three-window delta.
    flanks = ["".join(chars) for chars in itertools.product(ALPHABET, repeat=2)]
    for entry in rows:
        for left in flanks:
            for right in flanks:
                assert subtract(counts(left + entry["mutant"] + right),
                                counts(left + entry["reference"] + right)) == entry["delta"]
                outside_checks += 1
    # Literal fixed-length 200nt full-feature checks, covering all local edits.
    fixed_length_checks = 0
    for entry in rows:
        reference = "A" * 97 + entry["reference"] + "C" * 98
        mutant = "A" * 97 + entry["mutant"] + "C" * 98
        assert len(reference) == len(mutant) == 200
        assert subtract(counts(mutant), counts(reference)) == entry["delta"]
        fixed_length_checks += 1
    normalized_norms = {vector: sum(Fraction(abs(x), 198) for x in vector[:32])
                        + Fraction(abs(vector[-1]), 200) for vector in classes}
    smallest_norm = min(value for vector, value in normalized_norms.items() if any(vector))
    args.out.mkdir(parents=True)
    columns = ["delta_count_" + word for word in CANONICAL] + ["delta_gc_count"]
    write_csv(args.out / "all-snv-deltas.csv", ["reference_5mer", "alternative_center",
        "mutant_5mer", "rc_pair", "class_signature", *columns],
        ([e["reference"], e["alternative"], e["mutant"], e["rc_pair"],
          signature(e["delta"]), *e["delta"]] for e in rows))
    write_csv(args.out / "delta-classes.csv", ["class_signature", "oriented_edit_count",
        "rc_pair_count", "trimer_count_l1", "gc_count_l1", "normalized_feature_l1", *columns],
        ([signature(vector), len(entries), len({e["rc_pair"] for e in entries}),
          sum(map(abs, vector[:32])), abs(vector[-1]), str(normalized_norms[vector]), *vector]
         for vector, entries in sorted(classes.items())))
    write_csv(args.out / "zero-delta-edits.csv", ["reference_5mer", "alternative_center",
        "mutant_5mer", "rc_pair"],
        ([e["reference"], e["alternative"], e["mutant"], e["rc_pair"]] for e in zeros))
    pairs = [pair_checks(distance) for distance in (1, 2, 3)]
    summary = {"independent_implementation": True,
        "experimental_outcomes_read": False, "models_fitted": False,
        "arithmetic": "exact signed integer counts; Fraction for normalized L1",
        "canonical_trimers": CANONICAL, "sequence_length": 200,
        "trimer_frequency_denominator": 198, "gc_frequency_denominator": 200,
        "interior_positions_zero_based": [2, 197],
        "oriented_edits": len(rows), "reverse_complement_pairs": len(rc_pairs),
        "unique_delta_vectors": len(classes), "all_zero_edits": len(zeros),
        "all_zero_rc_pairs": len({e["rc_pair"] for e in zeros}),
        "zero_iff_mutated_5mer_equals_reference_reverse_complement": True,
        "class_size_histogram": dict(sorted(collections.Counter(map(len, classes.values())).items())),
        "trimer_count_l1_histogram_by_edit": dict(sorted(collections.Counter(
            sum(map(abs, e["delta"][:32])) for e in rows).items())),
        "minimum_nonzero_trimer_count_l1": min(sum(map(abs, v[:32])) for v in classes if any(v)),
        "minimum_nonzero_joint_integer_l1": min(sum(map(abs, v)) for v in classes if any(v)),
        "minimum_nonzero_normalized_feature_l1": str(smallest_norm),
        "minimum_normalized_l1_examples": [{"reference": e["reference"], "mutant": e["mutant"]}
            for e in rows if normalized_norms[e["delta"]] == smallest_norm][:8],
        "checks": {"all_oriented_delta_and_rc_checks": orientation_checks,
                   "all_two_base_flanks_full_count_checks": outside_checks,
                   "fixed_200nt_full_count_checks": fixed_length_checks},
        "double_edit_checks": pairs,
        "general_background_argument": "Any 3mer not containing the edited site is unchanged and cancels. For an interior SNV the only affected windows have starts p-2,p-1,p; their union is exactly the centered 5mer. Outside nucleotides therefore cannot alter this delta, for arbitrary background and boundary placement with p in [2,L-3].",
        "general_distance_argument": "A 3nt window cannot contain two positions separated by at least3nt. Each affected window then contains at most one edit, so its mixed second difference is zero. GC is an additive base count. Thus every pair separated by >=3nt has exactly zero interaction in this feature vector and every fixed linear predictor of it, including at sequence edges.",
        "limits": ["Interior 5mer enumeration does not apply to the first/last2positions without truncating the affected windows.",
            "Distinct vectors need not yield distinct predictions under particular coefficients; collisions or zero effects can increase.",
            "No quantitative minimum predictor effect follows without assumptions on coefficients; zero interaction here is a representation limitation, not a biological assertion.",
            "A fixed affine standardizer and intercept do not remove zero changes or feature-delta equivalences.",
            "Count norms and normalized-frequency norms are separate quantities."]}
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    script = Path(__file__).resolve()
    (args.out / "implementation.py").write_bytes(script.read_bytes())
    manifest = {p.name: digest(p) for p in sorted(args.out.iterdir()) if p.is_file()}
    (args.out / "manifest.json").write_text(json.dumps({"files": manifest}, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in ("oriented_edits", "reverse_complement_pairs",
        "unique_delta_vectors", "all_zero_edits", "minimum_nonzero_normalized_feature_l1",
        "checks", "double_edit_checks")}, indent=2))


if __name__ == "__main__":
    main()
