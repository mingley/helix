# Independent GC/trimer identity check

This check independently reconstructs the companion
[GC/trimer analysis](../gc-trimer-identity-v1/README.md), without importing its
implementation. It reads the same frozen sequence panel and counts GC bases
directly in each overlapping **raw** 3-mer, before any reverse-complement
canonicalization. The JSON parser also instantiates previously explored panel
observation fields, but the check accesses only each record's DNA sequence.
No new measurements or outcome values are used.

Every one of the **55,276** length-200 sequences satisfies the integer identity
`3 G = T + B`, where `G` counts all GC bases, `T` sums GC bases over the 198
raw 3-mer windows, and `B` is twice the GC count at the two ends plus the GC
count at their immediate neighbors. The check matches **all fields** of the
primary summary, including all seven boundary counts, all five correlation
sufficient statistics, the exact squared-correlation fraction, and its
30-digit decimal correlation `0.999575814973870283865674138435`. The decimal
is independently reconstructed from the square root of the exact fraction.

The retained 200-base examples have the same 198 canonical 3-mers, but GC
counts 66 and 67 and boundary values 0 and 3. Thus the original GC feature can
distinguish sequences with identical trimer features. These are constructed
sequences, not observed activity or biological effects.

## Exact coefficient relationship

Let a fixed sequence length be `L >= 4`, let `f_c` be canonical 3-mer class
`c`'s frequency, `k_c` its number of GC bases, `g = G/L`, and let `B` be the
boundary count above. Reverse complementation preserves each 3-mer's GC count,
so the raw-window identity implies

```text
g = ((L - 2) / (3L)) * sum_c(k_c * f_c) + B / (3L).
```

For arbitrary trimer coefficients `w_c` and GC coefficient `v`, the same
fixed score can therefore be expressed using `B` in place of GC by replacing
each `w_c` with `w_c + v * (L - 2) * k_c / (3L)` and giving `B` coefficient
`v/(3L)`. The verifier checks coefficient arithmetic for all 64 raw trimers
at lengths 4, 5 and 200, and independently evaluates nonuniform rational
scores for both constructed sequences and a panel example from each of the
seven boundary classes. The synthetic tests exercise all binary GC patterns of
lengths 4 through 8.

The boundary term matters: `0 <= B <= 6`, so omitting it changes the GC
fraction by up to `2/L`; it is not identically zero. Equal trimer features
alone do not force equal GC. The coefficient transformation preserves an
already specified score; it does not imply that refitting a ridge model with a
new feature basis gives identical fitted coefficients or predictions, because
the penalty changes. These arithmetic results make no claim of a new biological
mechanism or independent experimental replication.

## Reproduction

The [receipt](receipt.json) binds the panel, primary implementation, primary
summary and independent verifier by SHA-256. From the repository root:

```bash
python3 -m unittest discover -s research/variant-locality/gc-trimer-independent-v1 -p 'test_*.py'
gc_review_dir=$(mktemp -d)
python3 research/variant-locality/gc-trimer-independent-v1/verify.py \
  --receipt "$gc_review_dir/receipt.json"
cmp research/variant-locality/gc-trimer-independent-v1/receipt.json "$gc_review_dir/receipt.json"
```

The receipt path must be new; exclusive file creation prevents overwriting an
input or previous receipt. Synthetic controls check boundary extremes, invalid
DNA, the witness, exact rational reparameterization, an end-to-end altered
summary, and output-path safety under Python optimization. The two implementations share
the frozen panel and the stated identity; agreement is computational verification
of this analysis, not independent biological validation.
