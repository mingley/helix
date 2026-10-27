# Exact DNA feature-equivalence scan on the fixed panel

This outcome-free scan enumerates all **32,502,288** interior single-base substitution opportunities across the unchanged **55,276** deposited 200nt panel records. **745,052 opportunities (2.2923063139431907%)** leave all 32 reverse-complement-canonical overlapping 3mer counts and the GC numerator exactly unchanged. The representation therefore cannot distinguish the two sequences. These are proposed edits, not experimentally measured variant effects or evidence of biological neutrality.

All 3,072 directed reference-5mer/alternate classes occur in this panel and produce 1,405 distinct exact 33-integer delta vectors. Exactly 64 classes have an all-zero vector: the alternate is the complement of the reference center, and the four flanking bases are reverse-complement symmetric. Equivalently, the mutant 5mer is the reverse complement of the reference 5mer. Normalizing counts by 198 and GC by 200 preserves exact feature equality for these same-length sequences.

The number of zero-delta opportunities per panel record has median 13, mean 13.478761125985962, interquartile range 11–16, range 0–33, and descriptive type-7 2.5/97.5 percentiles 6–22. At least one occurs in 55,270 records; six have none. These are fixed-panel summaries, not inferential intervals.

`opportunity-classes.csv` retains every directed local class, its 33 exact deltas, panel-slot multiplicity and an actual representative ID/position. `per-element.csv` retains all record counts. `zero-feature-examples.csv` gives one actual representative for each of the 64 zero-delta classes, with source-sequence hash and both zero- and one-based core positions. Positions 2–197 zero-based (3–198 one-based) are included, giving 196 positions times three alternate bases per record.

The SHA-bound source is `datasets/context-panel/panel.json`. Only record IDs and DNA sequences are selected or used. Reporter responses, source effects, fitted models, clinical labels and statistical significance do not enter this computation. Every original record is retained: slots are not deduplicated or described as genomic-unique variants. There are 55,276 distinct exact DNA strings and 55,275 reverse-complement DNA identities in the input.

`checks.json` records exact reverse-complement and allele-swap identities for all 3,072 local classes, the exhaustive 64-class characterization, and direct full-200nt mutation feature recomputation for one actual representative of every observed class. Every exact check passes. `manifest.json` binds the six computational artifacts, including the exact `scan.py` implementation. `replay-verification.json` records a successful byte-for-byte replay of those six artifacts and their manifest; the receipt and this explanatory README are outside that six-file computation manifest.

Replay from the repository root into a new directory (Python standard library only):

```bash
python3 research/variant-locality/panel-scan-v1/scan.py \
  --panel datasets/context-panel/panel.json \
  --out /tmp/helix-feature-panel-replay
```

The script refuses an existing completed output. All six computational artifacts are deterministic and replay byte-for-byte. This scan excludes boundary substitutions, indels and multiple edits. Predictors with positional, longer-context, strand-specific or additional features need not share this limitation. The exact feature equivalence is a representation property; its quantification on this fixed panel is a bounded computational result, not a new biological mechanism.
