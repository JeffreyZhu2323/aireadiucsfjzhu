# AI-READI — surveys, CGM, wearables and cognition in a dysglycemia cohort

UCSF analysis workspace for the AI-READI dataset. All statistical testing is
complete; this README is the handoff document for writing the paper. Every
number below is reproduced from a committed result file, and each section names
the file and the script that produced it. `docs/RESULTS_INDEX.md` indexes every
result file individually if you need to locate a specific number.

**Cohort:** 2,245 participants, age 40–94 (mean 60.8, SD 11.2), in four
pre-assigned study groups.

| Study group (`meta_study_group`) | n | Median mean-CGM (mg/dL) |
|---|---:|---:|
| `healthy` | 761 | 115.6 |
| `pre_diabetes_lifestyle_controlled` | 552 | 122.4 |
| `oral_medication_and_or_non_insulin_injectable_medication_controlled` | 681 | 139.3 |
| `insulin_dependent` | 251 | 168.6 |

**Coverage of the primary measures:** HbA1c 2,182/2,245 · mean CGM 2,245/2,245 ·
MoCA 2,234/2,245 · median 11 days of CGM wear · 8 participants flagged
`cgm_quality_suspect`.

---

## The five findings

### 1. Psychosocial and socioeconomic burden tracks glycemia, but most of that signal is diabetes-severity confounding

This is the central result and the one that most shapes how the paper should be
framed.

Spearman correlations of 44 survey features against each glycemic outcome, BCa
bootstrap 95% CI (10,000 resamples, seed 42), Benjamini-Hochberg FDR applied
separately within each outcome. **26 of 44 features are FDR-significant against
HbA1c, and 26 of 44 against mean CGM.**

Strongest associations against HbA1c
(`outputs/survey_glycemia/survey_hba1c_spearman_bootstrap.csv`):

| Feature | rho | 95% CI | q | n |
|---|---:|---|---:|---:|
| `paid5_total` (diabetes distress) | +0.252 | [+0.211, +0.292] | 1.8e-30 | 2127 |
| `dm_selfcare_healthcare` | +0.236 | [+0.198, +0.275] | 2.9e-26 | 2065 |
| `sdoh_education_years` | −0.162 | [−0.204, −0.118] | 3.7e-12 | 1971 |
| `dm_selfcare_footcare` | +0.161 | [+0.121, +0.202] | 5.5e-13 | 2176 |
| `dm_selfcare_exercise` | −0.155 | [−0.197, −0.113] | 3.7e-12 | 2168 |
| `food_insecurity_score` | +0.155 | [+0.111, +0.197] | 3.7e-12 | 2146 |
| `sdoh_college_graduate` | −0.158 | [−0.202, −0.114] | 9.0e-12 | 1971 |
| `diet_sugary_drinks` | +0.117 | [+0.074, +0.159] | 2.8e-07 | 2094 |
| `sdoh_adversity_index` | +0.114 | [+0.069, +0.159] | 1.3e-06 | 1967 |
| `cesd10_depressed` | +0.072 | [+0.028, +0.115] | 1.5e-03 | 2179 |

Then the correction. `survey_glycemic_severity_stratified_pooling.py` re-estimates
every crude-significant association within each of the four study groups and
pools by inverse variance
(`outputs/survey_glycemia/survey_glycemic_severity_stratified_pooling.csv`):

- **Median crude rho falls 62% after stratification.** 49 of 52 feature-outcome
  pairs meet the ≥10% confounding threshold.
- Only **14 of 26** HbA1c associations and **7 of 26** mean-CGM associations
  survive pooled FDR.
- Roughly 76% of the median crude rank covariance is *between* severity strata
  rather than within them.

Survivors, with the magnitude of attenuation:

| Feature (HbA1c) | crude rho | pooled rho | change | pooled q |
|---|---:|---:|---:|---:|
| `paid5_total` | +0.252 | +0.099 | −61% | 1.4e-04 |
| `diet_sugary_drinks` | +0.117 | +0.087 | −26% | 8.9e-04 |
| `diet_risk_score` | +0.100 | +0.085 | −15% | 9.0e-04 |
| `dm_selfcare_exercise` | −0.155 | −0.070 | −55% | 4.8e-03 |
| `sdoh_education_years` | −0.162 | −0.066 | −59% | 1.1e-02 |
| `sdoh_college_graduate` | −0.158 | −0.067 | −58% | 1.1e-02 |
| `food_insecurity_score` | +0.155 | +0.062 | −60% | 1.1e-02 |
| `food_insecure` | +0.140 | +0.053 | −62% | 2.9e-02 |
| `dm_selfcare_footcare` | +0.161 | +0.052 | −68% | 3.1e-02 |

Casualties — crude-significant, pooled-null — include `dm_selfcare_healthcare`
(+0.236 → +0.003, −99%, q = 0.94), `sdoh_rx_cost_barrier_self` (−99%),
`sdoh_hospitalized_overnight` (−91%), and `cesd10_depressed` (−70%, q = 0.39).

**Two interpretive points the write-up should not miss.** First, the diabetes
self-care scales are coded *higher = better*, so the positive crude correlations
for `dm_selfcare_healthcare` and `dm_selfcare_footcare` mean that participants
reporting **more** healthcare engagement and footcare had **higher** HbA1c. That
is reverse causation operating through severity — sicker patients are told to
check their feet and attend clinic — and it is exactly the association that
collapses to zero under stratification. `dm_selfcare_exercise` behaves in the
expected direction and does survive. Second, `dm_selfcare_healthcare` and
`sdoh_college_graduate` fail the Cochran Q homogeneity test, so their pooled
estimates should not be reported as single cohort effects.

### 2. HbA1c shows consistently stronger survey correlations than mean CGM, but the margin is small

Tested at two levels, because they answer different claims.

**Aggregate, with the survey instrument as the unit of analysis**
(`outputs/survey_glycemia/hba1c_vs_cgm_magnitude_comparison.csv`, console
statistics preserved in `reports/hba1c_vs_cgm_magnitude_comparison.txt`): mean
|rho| is 0.0805 for HbA1c versus 0.0635 for mean CGM, and the correlation is
weaker under CGM for 32 of 44 features. Five of six instruments point the same
way (PAID-5 is the exception at parity, ratio 1.004). Mean difference −0.0140 on
the raw |rho| scale; exact cluster sign-flip permutation p = 0.0625, which is at
the floor of what six clusters can produce; cluster bootstrap 95% CI
[−0.0237, −0.0055]; mixed model with a random intercept per instrument
−0.0164 (SE 0.0050, p = 0.0011). On the Fisher-z scale the permutation reaches
p = 0.0312, the hard floor for k = 6.

The honest reading, which the script itself states: with six clusters this is a
statement about **consistency of sign across instruments**, not about effect
size.

**Per feature, on the same participants**
(`outputs/survey_glycemia/hba1c_vs_cgm_per_feature_steiger.csv`): HbA1c and mean
CGM correlate with each other at rho ≈ 0.68, so these are overlapping dependent
correlations and the test is Williams' t (Steiger 1980, eq. 14), with a paired
bootstrap CI on the difference in |rho| as the distribution-free check. HbA1c has
the larger |rho| for 35 of 44 features, but **only 4 reach FDR significance** —
and three of those are the same construct:

| Feature | rho HbA1c | rho CGM | Δ\|rho\| | Williams t | q |
|---|---:|---:|---:|---:|---:|
| `sdoh_education_code` | −0.161 | −0.085 | +0.076 | −4.27 | 5.0e-04 |
| `sdoh_education_years` | −0.162 | −0.086 | +0.075 | −4.25 | 5.0e-04 |
| `sdoh_college_graduate` | −0.158 | −0.086 | +0.072 | −4.05 | 7.8e-04 |
| `sdoh_english_first_language` | +0.019 | −0.039 | — | +3.25 | 1.3e-02 |

The fourth is not an attenuation result: the two correlations have opposite
signs, and its bootstrap CI [−0.069, +0.049] includes zero. Williams' t and the
bootstrap disagree on 7 features overall; the script flags these and the
bootstrap is the one to trust, since Williams' t is derived for Pearson
correlations and is only an approximation for Spearman.

So the defensible claim is narrow: **educational attainment is the one construct
that measurably attenuates when the glycemic anchor is switched from HbA1c to
mean CGM.** A general "CGM attenuates survey associations" claim is not
supported per-feature.

### 3. Richer CGM phenotyping does not rescue the null survey features

A reasonable hypothesis is that features unrelated to *average* glucose might
still track glycemic *shape* — time in range, excursion frequency, sustained
episodes. `null_survey_cgm_feature_scan.py` tests this directly: the 16 survey
features null against both HbA1c and mean CGM, crossed with 22 wear-normalized
CGM features, 352 tests with global FDR
(`outputs/survey_glycemia/null_survey_cgm_feature_scan.csv`).

**7 of 352 survive, and all 7 belong to a single feature**,
`sdoh_english_first_language`, which is nominal and whose sign is not
interpretable. They point consistently at hypoglycemia exposure (moderate-hypo
sustained windows rho = +0.106, CI [+0.065, +0.145], q = 6.2e-04, n = 2020; time
below range rho = +0.102; spikes per day rho = −0.077). On stratified re-check
the association holds in healthy, oral-medication and insulin groups but is null
in pre-diabetes.

This is a clean negative result and belongs in the paper as such: depression
severity, medical discrimination, housing problems and the rest do not track CGM
shape either.

### 4. Diabetes severity tracks physiology, not volitional activity

Kruskal-Wallis across the four study groups with epsilon-squared effect sizes,
BH-FDR within metric family, Dunn post-hoc where the omnibus survives
(`outputs/group_differences/`).

Graded monotonically with severity, with large effects for glucose and small but
consistent effects for autonomic physiology:

| Metric | epsilon² | q | healthy → insulin |
|---|---:|---:|---|
| `avg_spike_duration_minutes` | 0.352 | 1.5e-169 | 36 → 283 min |
| `mean_cgm_overall_mg_dl` | 0.297 | 4.5e-144 | 115.6 → 168.6 mg/dL |
| `avg_cgm_per_spike_mg_dl` | 0.259 | 9.2e-125 | 151 → 169 mg/dL |
| `mean_stress_level` | 0.055 | 2.1e-22 | 47.1 → 63.1 |
| `mean_heartrate_bpm` | 0.051 | 5.4e-21 | 75.9 → 82.2 bpm |
| `mean_respiratory_rate_bpm` | 0.028 | 1.8e-11 | 14.51 → 14.95 bpm |
| `mean_sleep_hours` | 0.010 | 2.5e-04 | 5.53 → 5.10 h |

Null across groups: **daily steps** (eps² = 0.0015, q = 0.40), **active
calories** (0.0014, q = 0.40), **SpO2** (0.0034, q = 0.18) and **daytime
movement** (0.0033, q = 0.10). Medians are essentially flat — steps range only
5,244–5,776 across all four groups.

The contrast is the finding: heart rate, stress, respiration and sleep separate
the groups cleanly while *how much people actually move* does not.

One non-monotonic result needs careful wording. Spike **frequency** is highest in
the pre-diabetes and oral-medication groups and **lowest in the insulin-dependent
group** (2.82/day versus 3.91 for oral medication, Dunn z = 7.13,
q = 6.0e-12), because insulin-dependent participants sit above the 140 mg/dL
threshold for long continuous stretches that count as one long episode rather
than many. Their episodes are ~8× longer. Spike count and spike duration must be
interpreted together or the direction looks paradoxical.

### 5. Cognition declines with glycemic burden, but the association is age- and severity-dependent

MoCA against glycemic and wearable measures, Spearman with BCa bootstrap CIs
(10,000 resamples, seed 42), age strata 40–54 / 55–69 / 70+
(`outputs/cognition/`).

Cohort-wide:

| Exposure | n | rho | 95% CI | p |
|---|---:|---:|---|---:|
| `avg_spike_duration_minutes` | 2227 | −0.185 | [−0.225, −0.143] | 1.2e-18 |
| mean glucose | 2233 | −0.174 | [−0.212, −0.133] | 1.5e-16 |
| `avg_cgm_per_spike_mg_dl` | 2227 | −0.150 | [−0.191, −0.108] | 1.1e-12 |
| `mean_sleep_hours` | 2098 | +0.138 | [+0.094, +0.181] | 2.0e-10 |
| `mean_stress_level` | 1988 | −0.088 | [−0.132, −0.043] | 8.1e-05 |
| `mean_respiratory_rate_bpm` | 1987 | −0.080 | [−0.122, −0.036] | 3.8e-04 |
| `mean_heartrate_bpm` | 1991 | −0.075 | [−0.119, −0.030] | 7.9e-04 |

Adjusted for age and study group (OLS, HC3 robust SEs, n = 2227, R² = 0.078),
spike burden remains independently associated with MoCA: `avg_spikes_per_day`
−0.134 (p = 1.3e-04) and `avg_cgm_per_spike_mg_dl` −0.0148 (p = 0.0089), while
`avg_spike_duration_minutes` drops to null (p = 0.680). Insulin-dependent status
carries −1.59 MoCA points (p = 1.9e-07) and oral medication −0.46 (p = 0.0053)
relative to healthy; age costs −0.049 points per year (p = 2.6e-16).

**The strata matter.** The glucose-cognition association is present at 40–54
(rho = −0.178) and 55–69 (−0.186) but **absent at 70+** (−0.060, p = 0.166), and
absent in the **insulin-dependent** group (rho = +0.0007, p = 0.991). The same
pattern repeats for spike features and every wearable stream. None of the 12
diabetes × age cells survives FDR.

Physical activity is largely an age artifact: the cohort-wide movement
correlations (total movement rho = +0.075) fall to rho ≈ 0.02 and p > 0.24 after
partialling out age. The one measure that survives age adjustment is
`physical_activity_num_days`, the number of days the monitor was worn
(partial rho = +0.061, p = 0.0049) — which is plausibly an adherence proxy rather
than an activity measure, and should be framed that way. Participants below the
MoCA < 26 cutoff did show lower movement and fewer wear days, but effect sizes
are small (rank-biserial −0.06 to −0.08).

Note that MoCA is **cross-sectional**: 2,254 participants have a MoCA score and
**none has a repeat test**, so no within-person cognitive change can be analyzed.

---

## Repository map

```
src/                        37 analysis scripts, flat, grouped by filename prefix
  add_*                     build the master table (run these first, in order)
  survey_*, hba1c_vs_cgm_*  finding 1 and 2
  null_survey_cgm_*         finding 3
  diabetes_group_*          finding 4
  moca_*, test2_, test3_    finding 5
  glucose_*, spike_*        CGM spike pipeline and its validation
  wearable_*, plot_sleep_*  wearable aggregation and sleep figures

outputs/
  cohort/                   master analysis table + survey codebook
  survey_glycemia/          findings 1-3, plus reports/ for console-only statistics
  group_differences/        finding 4 (omnibus + Dunn post-hoc pairs)
  cognition/                finding 5
  glycemic_phenotype/       spike episodes, person aggregates, method comparison
  figures/                  sleep interval and overlap charts
  eda_advanced/             partner-supplied outputs (see provenance warning below)
  archive/exploratory_2026-05-01/   superseded exploratory work, kept for provenance

docs/
  RESULTS_INDEX.md          every result file: shape, question, script, key columns
  survey_features.md        survey scoring rules, item mappings, data-quality fixes
  cgm_range_features.md     glycemic range definitions, sentinel handling, validation
  survey_instruments/       the five source survey PDFs
  ucsf_aireadi_project_notes.pdf

data/
  manifests/                AI-READI file manifests (committed)
  raw/omop/                 the five OMOP tables (NOT committed - see below)
```

The master table `outputs/cohort/full_person_aggregates.csv` is 2,245 × 97 and
is what most analyses read: one row per `person_id`, carrying spike metrics,
overall CGM, 25 CGM range features, 5 CGM quality flags, HbA1c, 10 wearable
streams, MoCA, 44 survey features, and study group and age.

## Getting the data and reproducing

The AI-READI source data is **not in this repository** and cannot be — it is
controlled-access, and `observation.csv` alone is 109 MB, above GitHub's 100 MB
per-file limit. Download the five OMOP tables from FAIRhub into
`data/raw/omop/` (`condition_occurrence`, `measurement`, `observation`,
`person`, `visit_occurrence`), and place the partner CGM export at
`outputs/glycemic_phenotype/glucose_raw.csv` (1.4 GB, ~6.2 M readings).

Scripts resolve paths relative to the repository, so it can be cloned anywhere.
Run from the repository root with the project interpreter:

```bash
.venv/bin/python src/glucose_spike_detection.py                        # spike episodes + person aggregates
.venv/bin/python src/wearable_stream_aggregates.py                     # wearable means
.venv/bin/python src/merge_wearable_into_full_person_aggregates.py     # start the master table
.venv/bin/python src/add_mean_cgm_overall_to_full_person_aggregates.py
.venv/bin/python src/add_cgm_ranges_to_full_person_aggregates.py
.venv/bin/python src/add_hba1c_to_full_person_aggregates.py
.venv/bin/python src/add_surveys_to_full_person_aggregates.py
```

Then any analysis script independently. All bootstraps are seeded (42, except
the impairment-activity family at 12345), so results are reproducible.
`src/hba1c_vs_cgm_per_feature_steiger.py` is the slowest at roughly 7 minutes.

## Cautions for the write-up

**Statistical**

- Every correlation in findings 1–3 is **crude**. Finding 1 quantifies how much
  of it is severity confounding; do not quote a crude rho without that context.
- FDR scope differs by analysis — within outcome for the survey scans, within
  metric family for the group differences, across strata for the MoCA glucose
  bootstrap, across metrics within stratum elsewhere. Match the claim to the
  correction that was actually applied.
- With n ≈ 2,200, significance is cheap. Epsilon-squared and the CIs carry the
  real information; several FDR-significant movement effects have eps² ≈ 0.006.
- No age or sex adjustment in the group-difference tests.
- Williams' t is approximate for Spearman; where it disagrees with the paired
  bootstrap, prefer the bootstrap.
- Columns ending `_items_answered` are completeness counts, not constructs. Some
  reach significance and should be reported as QC, not findings.

**Measurement**

- HbA1c and mean CGM are different constructs over different time windows — a
  single lab integral versus a ~11-day wearable mean. The comparison in finding 2
  is descriptive.
- `total_spikes` and the movement totals are confounded by wear duration. Prefer
  `avg_spikes_per_day` and the means-only variants.
- CGM sentinel readings are mapped Low → 39 and High → 401 mg/dL, which
  compresses both tails. The severe-hyper means are censored; `docs/cgm_range_features.md`
  recommends filtering `cgm_high_sentinel_frac > 0.10` (33 participants) for
  those columns specifically.
- PAID-5 asks diabetes-specific questions, yet 736 of 761 healthy participants
  answered it. Interpret with group restriction or adjustment.
- The exported `dietscore` column was corrupt; use the recomputed
  `diet_risk_score`. CES-D reverse-scoring was already applied at capture.
- `sdoh_job_loss_likely` has low coverage (n ≈ 859) because code 0 means
  "don't know". Uninsured and no-usual-care flags are sparse at ~2%.
- Spike detection uses a 140 mg/dL threshold with a one-reading dip tolerance
  across gaps ≤ 10 min. Interpolating episode boundaries does **not** change
  detection — both methods find 85,993 episodes cohort-wide, validated over
  5,000 synthetic streams and 40 real participants with zero mismatches. Only the
  ≥10 min minimum-duration filter changes counts, removing 10.2% of episodes.
  The master table uses the point-timestamp pipeline.

**Provenance — needs a decision before submission**

Every MoCA analysis reads
`outputs/eda_advanced/tables/master_patient_features.csv`, which was produced by
a partner script (`compiler_eda_advanced.py`) that is **not in this repository**.
That table cannot currently be regenerated from source, and
`outputs/eda_advanced/reports/walkthrough.md` still contains absolute paths from
the machine that produced it. Finding 5 therefore rests on an input whose
derivation is undocumented here. Either obtain the script from the partner or
rebuild those features from the OMOP tables before the methods section is
finalized. Note also that the partner walkthrough uses age bands
(<40, 40–55, >55) that do **not** match the production analyses (40–54, 55–69,
70+), and the <40 band is empty since cohort minimum age is 40.
