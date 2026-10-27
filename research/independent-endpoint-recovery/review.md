# Independent endpoint-recovery review

Independent parsing and arithmetic are complete. Reproducing a negative audit
means reproducing its computations and failures; it does not turn mismatching
source endpoints into agreement. All 36 original routes were estimable and
failed source agreement. Their arrays reproduced with maximum absolute
cross-implementation difference 4.3520742565306136e-14. The separate six-route
post-hoc audit also retained all estimable routes and all prescribed S12
derivatives; every source-agreement test failed. No route, population, value,
model or tolerance changed to improve agreement.

The explicit published-S11-ESS to full-S12 sample-standardization hypothesis
passed, including the fixed 55,276 admitted subset. It does not recover the
upstream activity-to-ESS recipe, certify producer history, add an independent
biological cohort or authorize changing previous targets. Original comparisons
use 1e-8 absolute plus 1e-8 relative tolerance; the separate S12 mapping uses
maximum absolute tolerance 1e-12. The six-route audit remains explicitly
post-hoc and uses its own separately frozen protocol.

Final commitments:

- Original independent manifest: `006e87a0876a9a002788b4367a6a88655cda6eed2fe5d0b54c9c687cf899a1ab`, retained inside [independent-recovery.tar.gz](../results/context-endpoint-recovery-v1/independent-recovery.tar.gz) at `independent-python/manifest.json`.
- [Original final comparison](../results/context-endpoint-recovery-v1/comparison/comparison.json): 7,008 checks, zero failures, SHA256 `9bddcb8da47c27d8252157cf8894cf332f5087f6562e54e23eec6ede0dec095c`.
- [Six-route independent manifest](../results/context-endpoint-recovery-v1/addenda/context-ess-cohort-audit-v1-independent/manifest.json): `38633a8390b8ef362c7a9a3c09a53bd731af51b61f04b6a2a8fc08b3e7b38b44`.
- [Six-route final comparison](../results/context-endpoint-recovery-v1/addenda/context-ess-cohort-audit-v1-comparison-v2/comparison.json): 1,162 checks, zero failures, SHA256 `5e0a680e048320549df419ea50e1aeaaee4ddb185004bc358fe3ec26314ef5b0`.
- Frozen final verifier SHA256: `5c6ed489a8ec51de04f9b26e3564272b25d8993a0dba6370f65462a35cc43485`.

Final comparisons bind both primary and independent manifests and their frozen
protocol. Full arrays, source ledgers and implementation snapshots remain in
the [publication archive](../results/context-endpoint-recovery-v1/README.md).
All 28 synthetic controls passed. The final verifier is frozen for the clean
workflow replay. The original empirical run retains its earlier implementation
snapshot; the final comparator adds format adapters and the separate audit.

Retained failed comparison attempts:

1. [First original comparison](../results/context-endpoint-recovery-v1/history/context-endpoint-recovery-comparison/failure.json) initially expected plain hashes, while the primary manifest stores `{bytes, sha256}` objects. Primary bytes were intact; the initial corruption interpretation was corrected. No source repair occurred.
2. [Second original comparison](../results/context-endpoint-recovery-v1/history/context-endpoint-recovery-comparison-v2/failure.json) initially expected a list of post-normalization parameters instead of a context-keyed dictionary. Only the format adapter changed.
3. [First six-route comparison](../results/context-endpoint-recovery-v1/history/context-ess-cohort-audit-v1-comparison/comparison.json) retained 42 failures caused solely by reserved missing-category label spelling. The final receipt explicitly maps independent `__UNMAPPED__` to primary `__UNMAPPED_CATEGORY__` and checks exact IDs and category counts. This changes no cohort, value or transformation.

All three failed attempts are preserved. Their corrections alter format
interpretation only; all original negative scientific findings remain.
