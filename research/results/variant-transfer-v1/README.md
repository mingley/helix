# Retained v1 source-gate failure

Both original source readers stop on the same literal `NA` in cell AO42972,
K562 variant `chr12:72652619:A:G`. All 241 selected-cell records agree, including
240 ordinary numeric float bit patterns from 120 fully decoded paired variants.
The frozen policy accepted no undocumented missing token. **The source gate
failed; no model-accuracy comparison or bootstrap was computed.**

`primary-source/`, `independent-source/` and `gate-comparison/` preserve every
original artifact and manifest. The independent reader exported fixed predictions
before the source failure; those are retained and are not performance evidence.
Openpyxl may instantiate other worksheet values internally. The separately
versioned exploratory v2 does not change this v1 failure.

See the [study report](../../variant-transfer-study.md) for source/assembly,
control/exposure exclusions, the two committed freezes and all limitations.
The final metadata cohort had 2,172 variants in 172 original components.

Verify without new source access or dependencies:

```bash
python3 -S research/results/variant-transfer-v1/packaging-implementation.py \
  --verify research/results/variant-transfer-v1
```

`publication.json` inventories all published bytes except itself. It binds both
original run manifests, their failed source comparison and the frozen authority.
