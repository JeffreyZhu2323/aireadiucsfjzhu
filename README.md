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

| Feature | Instrument | Direction | rho | 95% CI | q | n |
|---|---|---|---|---|---|---|
| `paid5_total` | PAID-5 | higher = worse | +0.252 | [+0.211, +0.292] | 1.8e-30 | 2127 |
| `dm_selfcare_healthcare` | Self-care | higher = better | +0.236 | [+0.198, +0.275] | 2.9e-26 | 2065 |
| `dm_selfcare_footcare` | Self-care | higher = better | +0.161 | [+0.121, +0.202] | 5.5e-13 | 2176 |
| `dm_selfcare_exercise` | Self-care | higher = better | -0.155 | [-0.197, -0.113] | 3.7e-12 | 2168 |
| `sdoh_education_years` | SDOH | higher = better | -0.162 | [-0.204, -0.118] | 3.7e-12 | 1971 |
| `food_insecurity_score` | Food insecurity | higher = worse | +0.155 | [+0.111, +0.197] | 3.7e-12 | 2146 |
| `sdoh_education_code` | SDOH | higher = better | -0.161 | [-0.203, -0.117] | 3.7e-12 | 1971 |
| `sdoh_college_graduate` | SDOH | higher = better | -0.158 | [-0.202, -0.114] | 9.0e-12 | 1971 |
| `food_insecurity_hvs_positive` | Food insecurity | higher = worse | +0.149 | [+0.105, +0.193] | 2.5e-11 | 2132 |
| `food_insecure` | Food insecurity | higher = worse | +0.140 | [+0.094, +0.182] | 3.5e-10 | 2146 |
| `diet_sugary_drinks` | Dietary | higher = worse | +0.117 | [+0.074, +0.159] | 2.8e-07 | 2094 |
| `paid5_distress` | PAID-5 | higher = worse | +0.112 | [+0.070, +0.157] | 7.4e-07 | 2127 |
| `sdoh_adversity_index` | SDOH | higher = worse | +0.114 | [+0.069, +0.159] | 1.3e-06 | 1967 |
| `sdoh_rx_cost_barrier_count` | SDOH | higher = worse | +0.108 | [+0.060, +0.154] | 7.3e-06 | 1904 |
| `diet_risk_score` | Dietary | higher = worse | +0.100 | [+0.057, +0.142] | 1.3e-05 | 2097 |
| `sdoh_insurance_type_count` | SDOH | nominal | +0.103 | [+0.055, +0.150] | 3.4e-05 | 1793 |
| `sdoh_rx_cost_barrier_self` | SDOH | higher = worse | +0.092 | [+0.045, +0.137] | 1.2e-04 | 1967 |
| `sdoh_urgent_care_visits_code` | SDOH | higher = worse | +0.083 | [+0.038, +0.129] | 5.4e-04 | 1967 |
| `sdoh_uninsured` | SDOH | higher = worse | +0.087 | [+0.036, +0.135] | 5.4e-04 | 1793 |
| `food_insecurity_items_answered` | Food insecurity | QC | -0.077 | [-0.117, -0.037] | 6.5e-04 | 2182 |
| `dm_selfcare_items_answered` | Self-care | QC | +0.077 | [+0.037, +0.115] | 6.8e-04 | 2182 |
| `cesd10_depressed` | CES-D-10 | higher = worse | +0.072 | [+0.028, +0.115] | 0.001 | 2179 |
| `sdoh_housing_instability` | SDOH | higher = worse | +0.069 | [+0.024, +0.112] | 0.004 | 1967 |
| `diet_fat_processed` | Dietary | higher = worse | +0.064 | [+0.021, +0.107] | 0.006 | 2099 |
| `sdoh_hospitalized_overnight` | SDOH | higher = worse | +0.053 | [+0.004, +0.100] | 0.033 | 1971 |
| `sdoh_care_access_barrier_count` | SDOH | higher = worse | +0.051 | [+0.005, +0.098] | 0.041 | 1940 |

And against mean CGM
(`outputs/survey_glycemia/survey_mean_cgm_spearman_bootstrap.csv`):

| Feature | Instrument | Direction | rho | 95% CI | q | n |
|---|---|---|---|---|---|---|
| `paid5_total` | PAID-5 | higher = worse | +0.223 | [+0.182, +0.263] | 2.5e-24 | 2190 |
| `dm_selfcare_healthcare` | Self-care | higher = better | +0.197 | [+0.157, +0.236] | 1.3e-18 | 2125 |
| `dm_selfcare_exercise` | Self-care | higher = better | -0.171 | [-0.211, -0.129] | 6.4e-15 | 2231 |
| `dm_selfcare_footcare` | Self-care | higher = better | +0.139 | [+0.100, +0.179] | 4.2e-10 | 2239 |
| `food_insecurity_hvs_positive` | Food insecurity | higher = worse | +0.116 | [+0.072, +0.159] | 4.4e-07 | 2193 |
| `food_insecurity_score` | Food insecurity | higher = worse | +0.110 | [+0.065, +0.152] | 1.8e-06 | 2207 |
| `paid5_distress` | PAID-5 | higher = worse | +0.107 | [+0.066, +0.149] | 3.1e-06 | 2190 |
| `food_insecure` | Food insecurity | higher = worse | +0.091 | [+0.044, +0.135] | 1.1e-04 | 2207 |
| `sdoh_uninsured` | SDOH | higher = worse | +0.096 | [+0.046, +0.142] | 1.9e-04 | 1841 |
| `sdoh_rx_cost_barrier_count` | SDOH | higher = worse | +0.086 | [+0.039, +0.131] | 5.8e-04 | 1956 |
| `diet_risk_score` | Dietary | higher = worse | +0.081 | [+0.038, +0.124] | 6.1e-04 | 2159 |
| `sdoh_rx_cost_barrier_self` | SDOH | higher = worse | +0.083 | [+0.037, +0.129] | 6.7e-04 | 2022 |
| `sdoh_urgent_care_visits_code` | SDOH | higher = worse | +0.081 | [+0.036, +0.125] | 8.5e-04 | 2024 |
| `sdoh_education_years` | SDOH | higher = better | -0.081 | [-0.123, -0.036] | 8.5e-04 | 2027 |
| `sdoh_college_graduate` | SDOH | higher = better | -0.080 | [-0.123, -0.035] | 8.7e-04 | 2027 |
| `sdoh_education_code` | SDOH | higher = better | -0.080 | [-0.122, -0.035] | 8.7e-04 | 2027 |
| `diet_sugary_drinks` | Dietary | higher = worse | +0.073 | [+0.029, +0.116] | 0.002 | 2156 |
| `food_insecurity_items_answered` | Food insecurity | QC | -0.069 | [-0.110, -0.027] | 0.003 | 2245 |
| `sdoh_insurance_type_count` | SDOH | nominal | +0.075 | [+0.028, +0.120] | 0.003 | 1841 |
| `dm_selfcare_items_answered` | Self-care | QC | +0.066 | [+0.027, +0.105] | 0.004 | 2245 |
| `diet_items_answered` | Dietary | QC | -0.061 | [-0.102, -0.020] | 0.008 | 2245 |
| `sdoh_adversity_index` | SDOH | higher = worse | +0.062 | [+0.017, +0.107] | 0.011 | 2022 |
| `cesd10_depressed` | CES-D-10 | higher = worse | +0.057 | [+0.014, +0.100] | 0.013 | 2242 |
| `diet_fat_processed` | Dietary | higher = worse | +0.058 | [+0.015, +0.100] | 0.013 | 2161 |
| `sdoh_hospitalized_overnight` | SDOH | higher = worse | +0.056 | [+0.010, +0.102] | 0.021 | 2028 |
| `dm_selfcare_diet` | Self-care | higher = better | -0.049 | [-0.091, -0.007] | 0.035 | 2212 |

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

| Feature | Target | crude rho | pooled rho | change | pooled 95% CI | pooled q | homog. |
|---|---|---|---|---|---|---|---|
| `paid5_total` | HbA1c | +0.252 | +0.099 | -61% | [+0.056, +0.141] | 1.4e-04 | yes |
| `sdoh_education_years` | HbA1c | -0.162 | -0.066 | -59% | [-0.109, -0.021] | 0.011 | yes |
| `sdoh_education_code` | HbA1c | -0.161 | -0.064 | -60% | [-0.108, -0.020] | 0.011 | yes |
| `dm_selfcare_footcare` | HbA1c | +0.161 | +0.052 | -68% | [+0.009, +0.093] | 0.031 | yes |
| `sdoh_college_graduate` | HbA1c | -0.158 | -0.067 | -58% | [-0.111, -0.023] | 0.011 | no |
| `dm_selfcare_exercise` | HbA1c | -0.155 | -0.070 | -55% | [-0.112, -0.028] | 0.005 | yes |
| `food_insecurity_score` | HbA1c | +0.155 | +0.062 | -60% | [+0.020, +0.104] | 0.011 | yes |
| `food_insecurity_hvs_positive` | HbA1c | +0.149 | +0.055 | -63% | [+0.013, +0.098] | 0.025 | yes |
| `food_insecure` | HbA1c | +0.140 | +0.053 | -62% | [+0.011, +0.095] | 0.029 | yes |
| `diet_sugary_drinks` | HbA1c | +0.117 | +0.087 | -26% | [+0.044, +0.129] | 8.9e-04 | yes |
| `paid5_distress` | HbA1c | +0.112 | +0.054 | -52% | [+0.012, +0.097] | 0.027 | yes |
| `diet_risk_score` | HbA1c | +0.100 | +0.085 | -15% | [+0.042, +0.127] | 9.0e-04 | yes |
| `food_insecurity_items_answered` | HbA1c | -0.077 | -0.074 | -4% | [-0.116, -0.032] | 0.003 | yes |
| `diet_fat_processed` | HbA1c | +0.064 | +0.079 | +22% | [+0.036, +0.121] | 0.002 | yes |
| `paid5_total` | mean CGM | +0.223 | +0.090 | -60% | [+0.048, +0.132] | 3.2e-04 | yes |
| `dm_selfcare_exercise` | mean CGM | -0.171 | -0.093 | -46% | [-0.134, -0.051] | 3.2e-04 | yes |
| `paid5_distress` | mean CGM | +0.107 | +0.064 | -41% | [+0.022, +0.105] | 0.019 | yes |
| `food_insecurity_items_answered` | mean CGM | -0.069 | -0.058 | -15% | [-0.099, -0.017] | 0.025 | yes |
| `diet_items_answered` | mean CGM | -0.061 | -0.057 | -7% | [-0.098, -0.016] | 0.025 | yes |
| `diet_fat_processed` | mean CGM | +0.058 | +0.060 | +4% | [+0.018, +0.102] | 0.025 | yes |
| `dm_selfcare_diet` | mean CGM | -0.049 | -0.071 | +44% | [-0.112, -0.029] | 0.008 | yes |
All 31 crude-significant associations that do **not** survive pooling:

| Feature | Target | crude rho | pooled rho | change | pooled 95% CI | pooled q | homog. |
|---|---|---|---|---|---|---|---|
| `dm_selfcare_healthcare` | HbA1c | +0.236 | +0.003 | -99% | [-0.041, +0.046] | 0.941 | no |
| `sdoh_adversity_index` | HbA1c | +0.114 | +0.038 | -67% | [-0.006, +0.082] | 0.135 | yes |
| `sdoh_rx_cost_barrier_count` | HbA1c | +0.108 | +0.022 | -79% | [-0.023, +0.067] | 0.393 | yes |
| `sdoh_insurance_type_count` | HbA1c | +0.103 | +0.050 | -51% | [+0.004, +0.096] | 0.059 | yes |
| `sdoh_rx_cost_barrier_self` | HbA1c | +0.092 | +0.001 | -99% | [-0.043, +0.045] | 0.966 | yes |
| `sdoh_uninsured` | HbA1c | +0.087 | +0.028 | -68% | [-0.019, +0.074] | 0.314 | yes |
| `sdoh_urgent_care_visits_code` | HbA1c | +0.083 | +0.040 | -51% | [-0.004, +0.085] | 0.113 | yes |
| `dm_selfcare_items_answered` | HbA1c | +0.077 | -0.014 | -82% | [-0.056, +0.028] | 0.587 | yes |
| `cesd10_depressed` | HbA1c | +0.072 | +0.022 | -70% | [-0.020, +0.064] | 0.386 | yes |
| `sdoh_housing_instability` | HbA1c | +0.069 | +0.043 | -38% | [-0.002, +0.087] | 0.096 | yes |
| `sdoh_hospitalized_overnight` | HbA1c | +0.053 | -0.005 | -91% | [-0.049, +0.040] | 0.906 | yes |
| `sdoh_care_access_barrier_count` | HbA1c | +0.051 | +0.031 | -39% | [-0.014, +0.076] | 0.237 | yes |
| `dm_selfcare_healthcare` | mean CGM | +0.197 | +0.001 | -99% | [-0.041, +0.044] | 0.983 | no |
| `dm_selfcare_footcare` | mean CGM | +0.139 | +0.039 | -72% | [-0.002, +0.081] | 0.165 | yes |
| `food_insecurity_hvs_positive` | mean CGM | +0.116 | +0.023 | -80% | [-0.018, +0.065] | 0.470 | yes |
| `food_insecurity_score` | mean CGM | +0.110 | +0.014 | -87% | [-0.028, +0.056] | 0.739 | yes |
| `sdoh_uninsured` | mean CGM | +0.096 | +0.056 | -41% | [+0.011, +0.102] | 0.051 | yes |
| `food_insecure` | mean CGM | +0.091 | -0.002 | -97% | [-0.044, +0.039] | 0.983 | yes |
| `sdoh_rx_cost_barrier_count` | mean CGM | +0.086 | +0.005 | -94% | [-0.039, +0.049] | 0.934 | yes |
| `sdoh_rx_cost_barrier_self` | mean CGM | +0.083 | +0.000 | -99% | [-0.043, +0.044] | 0.985 | yes |
| `diet_risk_score` | mean CGM | +0.081 | +0.050 | -38% | [+0.008, +0.092] | 0.057 | yes |
| `sdoh_urgent_care_visits_code` | mean CGM | +0.081 | +0.026 | -67% | [-0.017, +0.070] | 0.467 | yes |
| `sdoh_education_years` | mean CGM | -0.081 | +0.024 | -71% | [-0.020, +0.067] | 0.470 | yes |
| `sdoh_college_graduate` | mean CGM | -0.080 | +0.022 | -72% | [-0.022, +0.066] | 0.491 | yes |
| `sdoh_education_code` | mean CGM | -0.080 | +0.026 | -68% | [-0.018, +0.069] | 0.467 | yes |
| `sdoh_insurance_type_count` | mean CGM | +0.075 | +0.027 | -64% | [-0.019, +0.073] | 0.467 | yes |
| `diet_sugary_drinks` | mean CGM | +0.073 | +0.031 | -58% | [-0.011, +0.073] | 0.355 | yes |
| `dm_selfcare_items_answered` | mean CGM | +0.066 | -0.008 | -88% | [-0.050, +0.033] | 0.863 | yes |
| `sdoh_adversity_index` | mean CGM | +0.062 | -0.013 | -80% | [-0.056, +0.031] | 0.753 | yes |
| `cesd10_depressed` | mean CGM | +0.057 | +0.005 | -92% | [-0.037, +0.046] | 0.934 | yes |
| `sdoh_hospitalized_overnight` | mean CGM | +0.056 | -0.012 | -78% | [-0.056, +0.031] | 0.753 | yes |
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
weaker under CGM for 32 of 44 features.

| Instrument | features | mean \|rho\| HbA1c | mean \|rho\| CGM | mean diff | CGM/HbA1c |
|---|---|---|---|---|---|
| PhenX food insecurity | 4 | 0.1302 | 0.0962 | -0.0340 | 0.739 |
| PhenX SDOH | 20 | 0.0713 | 0.0487 | -0.0226 | 0.683 |
| CES-D-10 depression | 6 | 0.0307 | 0.0173 | -0.0134 | 0.562 |
| Dietary questionnaire | 5 | 0.0727 | 0.0626 | -0.0100 | 0.862 |
| Diabetes self-care | 6 | 0.1131 | 0.1085 | -0.0045 | 0.960 |
| PAID-5 diabetes distress | 3 | 0.1223 | 0.1228 | +0.0005 | 1.004 | Five of six instruments point the same
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

| Feature | n | rho HbA1c | rho CGM | diff \|rho\| | Williams t | Williams q | bootstrap 95% CI | Williams sig. |
|---|---|---|---|---|---|---|---|---|
| `sdoh_education_code` | 1971 | -0.161 | -0.085 | +0.076 | -4.27 | 5.0e-04 | [+0.039, +0.112] | yes |
| `sdoh_education_years` | 1971 | -0.162 | -0.086 | +0.075 | -4.25 | 5.0e-04 | [+0.039, +0.111] | yes |
| `sdoh_college_graduate` | 1971 | -0.158 | -0.086 | +0.072 | -4.05 | 7.8e-04 | [+0.036, +0.108] | yes |
| `sdoh_adversity_index` | 1967 | +0.114 | +0.065 | +0.049 | +2.76 | 0.051 | [+0.015, +0.084] | **no** |
| `food_insecurity_score` | 2146 | +0.155 | +0.112 | +0.043 | +2.51 | 0.088 | [+0.008, +0.077] | **no** |
| `diet_sugary_drinks` | 2094 | +0.117 | +0.076 | +0.042 | +2.36 | 0.108 | [+0.007, +0.075] | **no** |
| `food_insecure` | 2146 | +0.140 | +0.100 | +0.040 | +2.34 | 0.108 | [+0.007, +0.073] | **no** |
| `dm_selfcare_healthcare` | 2065 | +0.236 | +0.197 | +0.039 | +2.28 | 0.110 | [+0.006, +0.072] | **no** |
| `food_insecurity_hvs_positive` | 2132 | +0.149 | +0.114 | +0.035 | +2.03 | 0.169 | [+0.000, +0.070] | **no** |
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
interpretable. All seven, against `sdoh_english_first_language`:

| CGM feature | CGM family | n | rho | 95% CI | q |
|---|---|---|---|---|---|
| `cgm_mod_hypo_sustained_windows_per_day` | sustained excursions | 2020 | +0.106 | [+0.065, +0.145] | 6.2e-04 |
| `cgm_time_below_range_frac` | time in range | 2020 | +0.102 | [+0.061, +0.143] | 7.0e-04 |
| `cgm_mod_hypo_frac_time` | time in range | 2020 | +0.101 | [+0.059, +0.141] | 7.0e-04 |
| `cgm_mod_hypo_windows_per_day` | excursion frequency | 2020 | +0.097 | [+0.055, +0.137] | 0.001 |
| `avg_spikes_per_day` | spike pipeline | 2020 | -0.077 | [-0.118, -0.035] | 0.038 |
| `cgm_sev_hypo_windows_per_day` | excursion frequency | 2020 | +0.074 | [+0.031, +0.113] | 0.048 |
| `cgm_sev_hypo_frac_time` | time in range | 2020 | +0.073 | [+0.031, +0.113] | 0.048 |
On stratified re-check the association holds in healthy, oral-medication and
insulin groups but is null in pre-diabetes.

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

Which groups actually differ, from Dunn post-hoc on the six pairwise contrasts
(`diabetes_group_differences_dunn.csv`, z statistics, **bold** = BH-significant
at 0.05):

| Metric | healthy vs pre | healthy vs oral | healthy vs insulin | pre vs oral | pre vs insulin | oral vs insulin |
|---|---|---|---|---|---|---|
| `avg_cgm_per_spike_mg_dl` | **-4.47** | **-16.80** | **-20.64** | **-11.11** | **-16.46** | **-8.36** |
| `avg_spike_duration_minutes` | **-7.29** | **-21.60** | **-23.18** | **-12.77** | **-16.82** | **-7.44** |
| `avg_spikes_per_day` | **-3.27** | **-3.84** | **+4.45** | -0.35 | **+6.65** | **+7.13** |
| `mean_heartrate_bpm` | -0.79 | **-7.26** | **-7.82** | **-5.93** | **-6.90** | **-2.49** |
| `mean_respiratory_rate_bpm` | +0.36 | **-4.00** | **-6.26** | **-4.04** | **-6.25** | **-3.29** |
| `mean_sleep_hours` | +2.14 | **+3.77** | **+3.75** | +1.40 | +2.02 | +0.99 |
| `mean_stress_level` | -0.77 | **-7.45** | **-8.18** | **-6.12** | **-7.26** | **-2.70** |
| `night_movement` | +0.53 | **+3.37** | **+2.82** | **+2.61** | **+2.32** | +0.36 |
| `total_movement` | +0.11 | **+2.52** | **+2.59** | **+2.23** | **+2.40** | +0.75 |
| `total_spikes` | **-3.43** | **-3.56** | **+4.82** | +0.06 | **+7.13** | **+7.30** |
Read the physiology rows left to right: healthy and pre-diabetes are generally
indistinguishable and the separation appears once medication begins. Sleep is
the weakest — only the two healthy-versus-medicated contrasts survive.

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

Adjusted for age and study group (OLS, HC3 robust SEs, n = 2227, R² = 0.078,
model p = 1.5e-33), spike burden stays independently associated with MoCA
(`moca_score_adjusted_regression_coefficients.csv`):

| Term | Coefficient | SE | z | p | 95% CI |
|---|---|---|---|---|---|
| Intercept | +31.8543 | 0.9341 | +34.10 | 7.1e-255 | [+30.0234, +33.6851] |
| insulin-dependent (vs healthy) | -1.5853 | 0.3042 | -5.21 | 1.9e-07 | [-2.1816, -0.9891] |
| oral medication (vs healthy) | -0.4609 | 0.1654 | -2.79 | 0.005 | [-0.7851, -0.1367] |
| pre-diabetes (vs healthy) | -0.1496 | 0.1620 | -0.92 | 0.356 | [-0.4670, +0.1679] |
| `avg_spikes_per_day` | -0.1339 | 0.0349 | -3.83 | 1.3e-04 | [-0.2023, -0.0654] |
| `avg_spike_duration_minutes` | +0.0000 | 0.0001 | +0.41 | 0.680 | [-0.0001, +0.0002] |
| `avg_cgm_per_spike_mg_dl` | -0.0148 | 0.0057 | -2.61 | 0.009 | [-0.0259, -0.0037] |
| `meta_age` | -0.0491 | 0.0060 | -8.19 | 2.6e-16 | [-0.0609, -0.0374] |
Spike *frequency* and in-spike glucose survive adjustment; spike *duration*
drops to null, so the strong univariate duration correlation is largely carried
by age and severity.

All five wearable streams are associated with MoCA cohort-wide
(`moca_wearable_stream_bootstrap.csv`):

| Feature | n | rho | 95% CI | q |
|---|---|---|---|---|
| `mean_heartrate_bpm` | 1991 | -0.075 | [-0.119, -0.030] | 9.8e-04 |
| `mean_oxygen_saturation_pct` | 1620 | +0.071 | [+0.023, +0.119] | 0.004 |
| `mean_stress_level` | 1988 | -0.088 | [-0.132, -0.043] | 2.0e-04 |
| `mean_respiratory_rate_bpm` | 1987 | -0.080 | [-0.122, -0.036] | 6.2e-04 |
| `mean_sleep_hours` | 2098 | +0.138 | [+0.094, +0.181] | 1.0e-09 |
**The strata matter, and this is the main qualifier on finding 5.** Mean glucose
versus MoCA, by age band and by study group:

| Stratum | Feature | n | rho | 95% CI | q | FDR sig. |
|---|---|---|---|---|---|---|
| 40-54 | `mean_glucose` | 713 | -0.178 | [-0.249, -0.104] | 2.4e-06 | yes |
| 55-69 | `mean_glucose` | 982 | -0.186 | [-0.246, -0.124] | 1.3e-08 | yes |
| 70+ | `mean_glucose` | 538 | -0.060 | [-0.143, +0.025] | 0.166 | **no** |
| Stratum | Feature | n | rho | 95% CI | q | FDR sig. |
|---|---|---|---|---|---|---|
| healthy | `mean_glucose` | 759 | -0.144 | [-0.214, -0.074] | 2.6e-04 | yes |
| pre-diabetes | `mean_glucose` | 549 | -0.081 | [-0.162, +0.003] | 0.078 | **no** |
| oral med | `mean_glucose` | 676 | -0.109 | [-0.184, -0.030] | 0.009 | yes |
| insulin | `mean_glucose` | 249 | +0.001 | [-0.126, +0.129] | 0.991 | **no** |
The association holds at 40–54 and 55–69 but is **absent at 70+**, and absent in
both the pre-diabetes and **insulin-dependent** groups. The same pattern repeats
for spike features (FDR-significant rows only):

| Stratum | Feature | n | rho | 95% CI | q |
|---|---|---|---|---|---|
| 40-54 | `avg_spike_duration_minutes` | 710 | -0.181 | [-0.252, -0.106] | 3.6e-06 |
| 40-54 | `avg_cgm_per_spike_mg_dl` | 710 | -0.134 | [-0.207, -0.061] | 5.2e-04 |
| 55-69 | `avg_spike_duration_minutes` | 980 | -0.204 | [-0.265, -0.140] | 3.4e-10 |
| 55-69 | `avg_cgm_per_spike_mg_dl` | 980 | -0.184 | [-0.244, -0.120] | 1.0e-08 |
| Stratum | Feature | n | rho | 95% CI | q |
|---|---|---|---|---|---|
| healthy | `avg_spike_duration_minutes` | 756 | -0.163 | [-0.232, -0.093] | 2.1e-05 |
| healthy | `avg_spikes_per_day` | 760 | -0.101 | [-0.173, -0.027] | 0.008 |
| healthy | `avg_cgm_per_spike_mg_dl` | 756 | -0.092 | [-0.160, -0.022] | 0.011 |
| oral med | `avg_spike_duration_minutes` | 676 | -0.110 | [-0.185, -0.036] | 0.013 |
| oral med | `avg_cgm_per_spike_mg_dl` | 676 | -0.094 | [-0.169, -0.019] | 0.022 |
| Stratum | Feature | n | rho | 95% CI | q |
|---|---|---|---|---|---|
| healthy / 40-54 | `avg_cgm_per_spike_mg_dl` | 273 | -0.145 | [-0.256, -0.031] | 0.049 |
| healthy / 55-69 | `avg_spike_duration_minutes` | 313 | -0.165 | [-0.276, -0.057] | 0.010 |
| healthy / 55-69 | `avg_spikes_per_day` | 314 | -0.123 | [-0.235, -0.007] | 0.044 |
And for the wearable streams:

| Stratum | Feature | n | rho | 95% CI | q |
|---|---|---|---|---|---|
| 40-54 | `mean_heartrate_bpm` | 647 | -0.149 | [-0.223, -0.072] | 3.6e-04 |
| 40-54 | `mean_stress_level` | 646 | -0.122 | [-0.194, -0.045] | 0.002 |
| 40-54 | `mean_respiratory_rate_bpm` | 646 | -0.125 | [-0.197, -0.047] | 0.002 |
| 40-54 | `mean_sleep_hours` | 677 | +0.146 | [+0.071, +0.221] | 3.6e-04 |
| 55-69 | `mean_heartrate_bpm` | 873 | -0.127 | [-0.193, -0.061] | 2.8e-04 |
| 55-69 | `mean_stress_level` | 871 | -0.141 | [-0.206, -0.074] | 7.6e-05 |
| 55-69 | `mean_sleep_hours` | 921 | +0.170 | [+0.105, +0.233] | 1.0e-06 |
| Stratum | Feature | n | rho | 95% CI | q |
|---|---|---|---|---|---|
| healthy | `mean_heartrate_bpm` | 686 | -0.098 | [-0.173, -0.019] | 0.017 |
| healthy | `mean_stress_level` | 686 | -0.107 | [-0.179, -0.031] | 0.013 |
| healthy | `mean_sleep_hours` | 717 | +0.154 | [+0.081, +0.227] | 1.8e-04 |
| pre-diabetes | `mean_oxygen_saturation_pct` | 373 | +0.126 | [+0.027, +0.224] | 0.038 |
| pre-diabetes | `mean_respiratory_rate_bpm` | 492 | -0.138 | [-0.223, -0.053] | 0.011 |
| oral med | `mean_sleep_hours` | 624 | +0.136 | [+0.057, +0.215] | 0.003 |
| Stratum | Feature | n | rho | 95% CI | q |
|---|---|---|---|---|---|
| healthy / 40-54 | `mean_heartrate_bpm` | 254 | -0.169 | [-0.285, -0.045] | 0.035 |
| healthy / 55-69 | `mean_heartrate_bpm` | 281 | -0.160 | [-0.279, -0.040] | 0.012 |
| healthy / 55-69 | `mean_stress_level` | 281 | -0.182 | [-0.292, -0.065] | 0.006 |
| healthy / 55-69 | `mean_sleep_hours` | 296 | +0.239 | [+0.129, +0.349] | 1.6e-04 |
| oral med / 40-54 | `mean_sleep_hours` | 162 | +0.212 | [+0.052, +0.355] | 0.034 |
Note what is *absent* from these tables: **no stratum of the insulin-dependent
group reaches significance for any exposure**, and none of the 12 diabetes × age
cells survives FDR for mean glucose. Every surviving cross-stratified cell is
either healthy or oral-medication.

**Physical activity is largely an age artifact.** The crude and age-partial
correlations side by side make this unambiguous
(`moca_physical_activity_correlation.csv`):

| Metric | n | crude rho | crude q | age-partial rho | age-partial p |
|---|---|---|---|---|---|
| `total_movement` | 2107 | +0.075 | 0.003 | +0.017 | 0.430 |
| `day_movement` | 2107 | +0.068 | 0.005 | +0.025 | 0.248 |
| `night_movement` | 2107 | +0.064 | 0.007 | +0.000 | 0.988 |
| `physical_activity_num_days` | 2130 | +0.060 | 0.009 | +0.061 | 0.005 |
| `average_daily_activity` | 2130 | +0.041 | 0.068 | -0.016 | 0.454 |
| `average_active_calories_kcal` | 2173 | +0.028 | 0.200 | -0.026 | 0.229 |
Three of the four crude-significant movement measures collapse to near zero once
age is partialled out. The only survivor is `physical_activity_num_days`, the
number of days the monitor was worn — plausibly an adherence or engagement proxy
rather than an activity measure, and it should be framed that way.

The impairment analysis points the same way
(`moca_impairment_activity_tests.csv`, FDR-significant rows):

| MoCA cutoff | Metric | n impaired / unimpaired | median impaired vs not | rank-biserial | 95% CI | q | mean-diff CI excl. 0 |
|---|---|---|---|---|---|---|---|
| 26 | `total_movement` | 853 / 1254 | 49,123 vs 52,114 | -0.075 | [-0.127, -0.023] | 0.010 | **no** |
| 26 | `day_movement` | 853 / 1254 | 26,895 vs 29,039 | -0.069 | [-0.120, -0.018] | 0.015 | **no** |
| 26 | `night_movement` | 853 / 1254 | 20,999 vs 22,475 | -0.062 | [-0.113, -0.013] | 0.022 | **no** |
| 26 | `physical_activity_num_days` | 860 / 1270 | 15 vs 16 | -0.082 | [-0.133, -0.033] | 0.007 | **no** |
| 24 | `physical_activity_num_days` | 459 / 1671 | 15 vs 15 | -0.081 | [-0.142, -0.019] | 0.045 | **no** |
| 23 | `physical_activity_num_days` | 329 / 1801 | 15 vs 15 | -0.092 | [-0.162, -0.023] | 0.044 | **no** |
Two cautions on this table. Effect sizes are small (rank-biserial −0.06 to
−0.09), and although the rank-biserial CIs exclude zero, **the mean-difference
CIs do not for any row** — this is a distributional shift, not a clean difference
in means. At cutoffs 24 and 23 the medians are identical (15 vs 15 days), so the
significance there is distribution shape rather than central tendency. Only
`physical_activity_num_days` survives the stricter cutoffs; the movement metrics
are significant at MoCA < 26 only.

Stratified, just two cells survive: `physical_activity_num_days` in the 55–69
band (rho = +0.092, q = 0.030, n = 930) and in 55–69 oral-medication
participants (rho = +0.219, q = 6.9e-04, n = 304).

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
