# Index of every result file

Companion to the findings in `README.md`. Use this to locate the file behind a
specific number. Shapes are rows × columns, excluding the header.

Shared conventions across result files: `spearman_rho` is Spearman's rho;
`ci_low`/`ci_high` are bootstrap 95% bounds; `p_value_fdr_bh` is the
Benjamini-Hochberg adjusted p-value and `significant_fdr_0_05` the verdict at
alpha = 0.05; `n` is the complete-case count for that specific test.

## cohort/ — the analysis tables everything else reads

| File | Shape | Contents |
|---|---|---|
| `full_person_aggregates.csv` | 2245 × 97 | **Master table.** One row per `person_id`: spike metrics, overall CGM, 25 CGM range features, 5 CGM quality flags, HbA1c, 10 wearable streams, MoCA, 44 survey features, study group, age. |
| `survey_codebook.csv` | 359 × 7 | Item × response-code map with labels and counts, for all six instruments. How to read any raw survey code. |
| `survey_items_by_person.csv` | 2280 × 97 | Item-level survey responses per person, before scale scoring. |
| `wearable_stream_means.csv` | 2184 × 11 | Per-person means of the wearable streams after sentinel cleaning. |

Built by `src/merge_wearable_into_full_person_aggregates.py` then the four
`src/add_*.py` scripts; surveys by `src/add_surveys_to_full_person_aggregates.py`
and `src/survey_codebook.py`. Scoring rules are in `docs/survey_features.md`;
CGM range definitions in `docs/cgm_range_features.md`.

## survey_glycemia/ — findings 1, 2 and 3

| File | Shape | Question it answers | Script |
|---|---|---|---|
| `survey_hba1c_spearman_bootstrap.csv` | 44 × 13 | Which survey features correlate with lab HbA1c? | `survey_glycemic_spearman_bootstrap.py` |
| `survey_mean_cgm_spearman_bootstrap.csv` | 44 × 13 | Same, against mean CGM. | same |
| `survey_glycemic_severity_stratified_pooling.csv` | 52 × 28 | **How much of each crude association is diabetes-severity confounding?** Carries per-stratum rho, inverse-variance pooled rho, Cochran Q, I², and the `survives` verdict. | `survey_glycemic_severity_stratified_pooling.py` |
| `hba1c_vs_cgm_magnitude_comparison.csv` | 6 × 7 | Across instruments, is HbA1c the stronger anchor? One row per instrument. | `hba1c_vs_cgm_magnitude_comparison.py` |
| `reports/hba1c_vs_cgm_magnitude_comparison.txt` | — | The permutation, cluster-bootstrap and mixed-model statistics for the above, which the script prints but does not write to CSV. | same |
| `hba1c_vs_cgm_per_feature_steiger.csv` | 44 × 21 | Per feature, is HbA1c stronger than mean CGM on the same participants? Williams' t plus a paired bootstrap on the difference in \|rho\|. | `hba1c_vs_cgm_per_feature_steiger.py` |
| `null_survey_cgm_feature_scan.csv` | 352 × 16 | Do the double-null survey features track CGM *shape* rather than average? 16 features × 22 CGM features. | `null_survey_cgm_feature_scan.py` |

Columns worth knowing: in the pooling file, `pct_change` is the shift from crude
to pooled rho and `confounded` flags a change of at least 10%; `homogeneous`
reports the Cochran Q test, and a false value means the pooled estimate should
not be quoted as a single cohort effect. In the Steiger file, `diff_abs` is
|rho_HbA1c| − |rho_CGM| (positive means HbA1c is stronger) while `diff_signed`
is the raw difference used by Williams' t; `signs_agree` false means the two
correlations point opposite ways and `diff_abs` is not an attenuation;
`boot_excludes_zero` is the distribution-free verdict and takes precedence where
it disagrees with `significant_fdr_0_05`.

## group_differences/ — finding 4

Kruskal-Wallis across the four study groups with epsilon-squared, BH-FDR within
metric family, and Dunn post-hoc on the six pairwise contrasts where the omnibus
survives. Each `*_dunn.csv` holds the post-hoc pairs for the omnibus file of the
same name.

| File | Shape | Metrics tested | Script |
|---|---|---|---|
| `diabetes_group_differences.csv` | 14 × 25 | 4 spike + 10 wearable, including exposure-weighted movement totals. | `diabetes_group_differences.py` |
| `diabetes_group_differences_means_only.csv` | 10 × 25 | 3 spike means + 7 wearable means; totals dropped to remove wear-duration confounding. **Prefer this one.** | `diabetes_group_differences_means_only.py` |
| `diabetes_group_differences_glucose_only.csv` | 4 × 25 | The 4 spike features alone. | `diabetes_group_differences_glucose_only.py` |
| `diabetes_group_differences_wearable_plus_mean_cgm.csv` | 8 × 25 | 7 wearable means + `mean_cgm_overall_mg_dl`. | `diabetes_group_differences_wearable_plus_mean_cgm.py` |

Each omnibus row carries median and IQR per group, `kruskal_H`, `epsilon_squared`
and `kruskal_p_fdr_bh`. Dunn rows carry `mean_rank_diff`, `z` and `p_value_bh`.

## cognition/ — finding 5

Age strata are 40–54 / 55–69 / 70+. Bootstraps are BCa, 10,000 resamples,
seed 42 (the impairment-activity family uses stratified percentile CIs, seed
12345). Impairment is MoCA below cutoff, evaluated at 26, 24 and 23.

| File | Shape | Question |
|---|---|---|
| `moca_glucose_bootstrap.csv` + `_by_age` / `_by_diabetes` / `_by_diabetes_age` | 1, 3, 4, 12 | Mean glucose vs MoCA, overall and per stratum. |
| `moca_spike_spearman_bootstrap.csv` + 3 stratified variants | 3, 9, 12, 36 | Spike aggregates vs MoCA. |
| `moca_wearable_stream_bootstrap.csv` + 3 stratified variants | 5, 15, 20, 60 | The five wearable vitals/sleep streams vs MoCA. |
| `moca_physical_activity_correlation.csv` | 6 × 12 | Activity vs MoCA cohort-wide, **including the age-partial Spearman** that shows most of it is an age artifact. |
| `moca_pa_correlation_by_age.csv` / `_by_diabetes` / `_by_age_diabetes` | 18, 24, 72 | Stratified activity correlations. |
| `moca_impairment_activity_tests.csv` | 18 × 33 | Impaired vs unimpaired activity: Mann-Whitney, rank-biserial effect size, bootstrap CIs. The definitive impairment test. |
| `moca_impairment_activity_by_*.csv` | 36, 48, 144 | Stratified Mann-Whitney versions. |
| `moca_impairment_activity_bootstrap_by_*.csv` | 54, 72, 216 | Stratified group-mean bootstrap CIs. No p-values — effect estimates only. |
| `moca_score_adjusted_regression_coefficients.csv` + `_model_summary.csv` | 8, 1 | OLS MoCA ~ spikes + age + study group, HC3 robust SEs. The multivariable result. |
| `moca_score_spearman_results.csv` | 7 × 6 | Legacy omnibus screen. Superseded by the dedicated bootstrap files above; keep only as a first-pass record. |

## glycemic_phenotype/ — the CGM spike pipeline

| File | Shape | Contents |
|---|---|---|
| `glucose_raw.csv` | 6,194,153 × 4 | Partner CGM export. **Not committed** (1.4 GB). Immutable input. |
| `glucose_spikes.csv` | 85,993 × 7 | One row per spike episode, point timestamps. The version behind the master table. |
| `glucose_spikes_interpolated.csv` | 85,993 × 8 | Same episodes with interpolated boundaries — identical count confirms interpolation does not change detection. |
| `glucose_spikes_interpolated_min10.csv` | 77,258 × 8 | After the ≥10 min filter: 10.2% fewer episodes. |
| `spike_aggregates_by_person.csv` (+ 2 variants) | 2245 × 6 | Per-person spike totals, rates, mean duration, mean in-spike glucose. |
| `spike_aggregate_method_comparison.csv` | 2245 × 10 | Per-person deltas between the point and min10 methods. |
| `spike_daily_by_cgm_day.csv` | 24,012 × 4 | Per person-day spike counts. |
| `glucose_summary_person_level.csv` / `_reading_level.csv` | 20 × 11 | Descriptive glucose summaries stratified by study group, at person and reading level. |

Spike definition, dip tolerance and the validation evidence are in `README.md`
under finding 4 and in `src/validate_spike_logic.py`.

## figures/

`sleep_interval_summary_chart.png` and `sleep_overlap_pairs_chart.png`. Both
plotting scripts read summaries that now live in
`outputs/archive/exploratory_2026-05-01/legacy_outputs/`, so regenerate from
there rather than from `outputs/` directly.

## Directories that are not paper results

- `eda_advanced/` — partner-supplied EDA. `tables/master_patient_features.csv`
  is a required input to every MoCA analysis but was produced outside this repo;
  see the provenance warning in `README.md`.
- `archive/exploratory_2026-05-01/` — exploratory work superseded on 2026-05-01,
  retained for provenance. `archive_actions.json` records what was moved and why.
- `manifest_lists/` — raw file path listings from `src/config.py`, a one-off
  utility.
