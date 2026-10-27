# Frozen regulatory study

This is a retrospective computational reanalysis of published reporter-effect estimates. It does not establish new biological mechanisms, independent laboratory replication, clinical risk, or whole-genome understanding.

## Fixed evaluation results

| Evaluation | Model | Macro position MSE | Macro position MAE |
| --- | --- | ---: | ---: |
| hbb_guarded_position_cv | mean | 0.021933886 | 0.105793703 |
| hbb_guarded_position_cv | chemistry | 0.022683072 | 0.109804500 |
| hbb_guarded_position_cv | coordinate | 0.026119801 | 0.122370989 |
| hbb_guarded_position_cv | delta3mer | 0.025027226 | 0.118152882 |
| hbb_guarded_position_cv | motif | 0.024996569 | 0.117568167 |
| hbb_guarded_position_cv | full | 0.029609654 | 0.129969335 |
| primary_leave_element_out | mean | 0.116392187 | 0.179463475 |
| primary_leave_element_out | chemistry | 0.116899210 | 0.180352538 |
| primary_leave_element_out | coordinate | 0.116528763 | 0.182018632 |
| primary_leave_element_out | delta3mer | 0.118370897 | 0.185745763 |
| primary_leave_element_out | motif | 0.116829487 | 0.181056745 |
| primary_leave_element_out | full | 0.117988875 | 0.187666960 |
| hbb_to_external_transfer | mean | 0.116723522 | 0.171381539 |
| hbb_to_external_transfer | chemistry | 0.117124800 | 0.174229599 |
| hbb_to_external_transfer | coordinate | 0.115579645 | 0.175647228 |
| hbb_to_external_transfer | delta3mer | 0.118532205 | 0.186719475 |
| hbb_to_external_transfer | motif | 0.116935604 | 0.174668577 |
| hbb_to_external_transfer | full | 0.117009723 | 0.187469569 |
| row_random_diagnostic | mean | 0.021462133 | 0.104526924 |
| row_random_diagnostic | chemistry | 0.021725836 | 0.106542351 |
| row_random_diagnostic | coordinate | 0.020792210 | 0.104677718 |
| row_random_diagnostic | delta3mer | 0.020941712 | 0.106586400 |
| row_random_diagnostic | motif | 0.021855343 | 0.107710758 |
| row_random_diagnostic | full | 0.020436205 | 0.104838825 |

Primary loss first averages alleles within position, positions within assay, then assays equally. HBB was already explored. Row-random results are a deliberately easier/leaky diagnostic, not primary validation. Input quality keeps every passing allele, including incomplete positions. No score p-value or count is a biological predictor.

## Same-sequence context diagnostics

| Relationship | Matched variants | Coefficient MSE | Conditional context-blind label MSE bound |
| --- | ---: | ---: | ---: |
| library_repeat | 943 | 0.033933404 | 0.008483351 |
| 48h_vs_24h_context | 1368 | 0.056943348 | 0.014235837 |
| GBM_vs_HEK293T_cell_context | 770 | 0.102880390 | 0.025720097 |
| GABPA_knockdown_context | 769 | 0.033854096 | 0.008463524 |
| scrambled_control_context | 769 | 0.027762029 | 0.006940507 |

These are finite source-estimate bounds, not biological irreducibility. Noise, library, time, cell and perturbation context can all contribute. A repeat in another element is not a gene-specific noise baseline. Overlapping TERT comparisons are not independent samples.

## Source forensics and validation

The source-anchored HBB audit maps 62 database equality descriptions to original deletion rows. They produce 53 distinct altered target sequences. Raw duplicate-label inconsistency disappears after the predefined tag filter: 35 qualified rows have 35 distinct sequences. Related representation issues have prior public reports; global novelty is unestablished.

Full predictions, fits, fold/guard identities and source-derived features are in study.json.gz. Bootstrap draws and spatial controls are in resampling.json. Their tails are conditional robustness diagnostics, not calibrated biological-significance probabilities. Ninety-nine controls and six-comparison Bonferroni adjustment cannot resolve a tail below 0.05.

Run the independent Python implementation and comparator before claiming numerical reproduction. Preserve failed attempts; exact-tie rank metrics can be unstable under floating reduction order even when continuous predictions/losses agree. Read the source literature/novelty ledger before labeling an observation new.
