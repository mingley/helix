# Exploratory variant-transfer continuation after the v1 source failure

The [v2 machine protocol](variant-transfer-v2.json) continues the same fixed-model
Siraj reporter transfer question under an explicitly new operational missingness
rule. **V1 remains a failed strict source import.** V2 is exploratory after
partial outcome exposure, and its repository freeze must precede further bulk
effect parsing. It is not outcome-naive validation or external preregistration.

Protocol SHA256:
`26861a28a1c049bcb6da66d17ff25bb13143f45cbf166a224af2a5aa27f2e04c`.
Parent v1 SHA256:
`811930c2d3a5ee763014ea611acda4d2c84404b7d5ae399c7843c62d0a828149`.
The machine classification is `exploratory_after_v1_source_failure`.

## What was exposed

The committed v1 schema allowed no undocumented missing token. The primary XML
reader recorded **240 selected numeric effect cells**, then stopped at the
ordinary shared-string **`NA` in AO42972**, for **chr12:72652619:A:G, K562**.
Its retained cell ledger contains 241 entries, including the blocking cell.
This was an intended source gate failure, not evidence of a successful v1
transfer test.

The independent v1 reader also stopped at that cell. Its XLSX iteration may have
internally decoded nonselected worksheet values; the total exposure cannot be
certified as only those 240 selected numeric cells. The independent attempt also
retained a fixed-prediction artifact. Thus the accurate statement is that no new
model was fitted and no observed-target performance comparison was established,
rather than claiming that no predictions had ever been produced.

Both complete failed attempts are copied byte-for-byte into the new binding's
[history](../research/variant-transfer/binding-v2/history/primary-v1/manifest.json).
The [operational evidence](../research/variant-transfer/binding-v2/operational-missingness-evidence.json)
pins their failures, source-cell ledgers and existing fixed-predictor output.
No new S4 effect value was read while preparing v2.

## The sole new operational rule

An ordinary source text cell resolving **exactly to `NA`**, without a formula,
is operationally **unreported with unknown reason**. Exclude a variant from
complete-case evaluation unless both K562 and HepG2 effects are ordinary finite
numeric cells. This is a new analyst-defined rule justified by the observed
token and disclosed as an assumption. It does not establish producer-defined
missingness semantics, a particular QC reason or missing-at-random behavior.

Preserve every affected cell, ID and original component. Do not impute a value,
substitute zero, infer biological inactivity, or filter by observed magnitude,
sign, activity, significance or model error. Retain the already exposed values
under the same complete-case rule; do not carve out a conveniently fresh subset.

Formulas, errors, booleans, absent/blank cells, nonfinite numbers, numeric
strings and every other text value remain fatal. Altered case or whitespace
does not turn another token into `NA`.

The binding exposes this policy as:

```json
{
  "allowed_tokens": ["NA"],
  "operational_missingness": {
    "interpretation": "literal_unreported_unknown_reason",
    "rule": "exclude_variant_unless_both_contexts_ordinary_finite_numeric",
    "source_semantics_established": false
  }
}
```

`operational_missingness` is a top-level object in the actual protocol and
binding; the example groups the relevant fields for reading.

## What remains fixed

The metadata cohort is the same **2,172 variants in 172 original components**.
Its compressed bytes are identical to v1, SHA256
`4d1fd2725ee0ed2952cf95cefa555ba4418c1338bdcd047a34d849d0b95ddc8b`.
No new source, sequence, coordinate conversion, technical-control rule,
previous-source guard or component relation is introduced. The independently
reproduced global grouping still passes **4,689 categories with zero failures**.

Complete-case exclusions never regroup the remaining variants or rename their
components. An original component with no complete variant supplies no loss
mean and is absent from the complete-case bootstrap population. Count those
empty components and all omitted members explicitly. The final sample size is
unknown at this freeze.

The following machine blocks are copied unchanged from v1: contexts, endpoint,
grouping, pretrained model, primary comparison, secondary diagnostics,
uncertainty, synthetic controls, independent reproduction and
`analysis_contract`. In particular:

- Use the same publisher S4, SHA256
  `723c4a0e1cfdb2c9c6ae55b3fa9206855dd95e0241319b9fc8d0831d3b42f9c2`.
  The endpoint remains the cell-specific DESeq2 Alt-minus-Ref by RNA-minus-DNA
  interaction in published delta log2 fold-change units. The exact selected
  columns remain AN for HepG2 and AO for K562.
- Use the same five frozen activity fits, artifact manifest SHA256
  `2ea476bc5413dc13486dec9f3d7e9162a50a0e47f1d83ce29efabc5f5cc5af9d`.
  The primary predictor averages the five fixed matched-cell I deltas; S/A
  slopes remain secondary. No refitting, calibration or sign choice is allowed.
- Compare I against zero using equal-component macro MSE, averaging the two
  contexts equally within a variant and variants equally within a component.
  Useful transfer still requires positive improvement and a positive lower
  endpoint of the same conditional interval.
- Retain 2,000 paired whole-component bootstrap draws, lexicographic component
  order, SplitMix64 seed `0x2026102801` and Type7 2.5%/97.5% quantiles.
  At least 100 complete variants and 10 original components remain required
  for an interval claim. All failed comparisons and negative results remain.

## Required next freeze and interpretation

The parent must commit this protocol and the
[new exact binding](../research/variant-transfer/binding-v2/binding.json) before
either reader resumes bulk effects. That binding is SHA256
`a54badc4645d7ac88c31b7722fbbd89aaefa2b995bc0a02befd15123d83d4fba`.
Only then may independent readers apply the exact operational rule, retain the
full type/missingness ledger and run the fixed comparison.

Results will concern this exploratory complete-case reporter panel. Missingness
can change that panel's composition. Numerical agreement between programs does
not validate the missingness assumption, establish endogenous causality or make
the analysis outcome-naive. Any further missing-token or analytical change
requires another explicitly exposed version; a v2 result cannot convert v1's
source gate into a pass.
