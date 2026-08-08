# Cognition & Diabetes Advanced EDA Walkthrough

I have fully implemented your requests into a new dedicated script: `compiler_eda_advanced.py`. The script was successfully executed, processing the entire array of 5-minute and hourly datasets. 

## Completed Objectives

- **MOCA Validation**: Evaluated `moca_score` across all patients. All valid patients had exactly one discrete MOCA test recorded, preserving temporal data integrity.
- **Bootstrapping Analysis**: Executed 2,000 resampled permutation instances per pairwise diabetes classification (e.g., Healthy vs Type 2) to calculate the mean differences in **MOCA Cognitive Score** and **Mean Blood Glucose**. 
  - Results output directly to `bootstrap_diff_moca_score.csv` and `bootstrap_diff_mean_glucose.csv`.
- **Summary Tables**: Processed the individual time-series CSVs from the `5min` and `hourly` omnibus subsets to create master participant summary tables (concatenated and aggregated by `participant_id`):
  - `all_patients_5min_summary.csv`
  - `all_patients_hourly_summary.csv`

## Generated Visualizations & Analysis

The following plots were generated outlining cognitive stratification and physical health profiles across patients. 

### Stratified Cognitive Assessment
We segmented the participants by both **Age Bracket** (`<40`, `40-55`, `>55`) and **Diabetes Diagnostic Type**. As requested, this illustrates how overall MoCA outcomes distribute within equivalent age ranges conditioned strictly on diabetes states.

![Cognition Stratified by Age and Diabetes Class](/Users/guardian/.gemini/antigravity/brain/a7afde04-87b0-4744-b59e-ea6493409fed/artifacts/stratified_cognition.png)

### Cognition vs Physical Activity & Sleep
We computed a Spearman correlation block mapping `moca_score`, overall daily activity metrics (`total_movement`, `day_hr`), and sleep metrics (`night_mean_hr`, `night_movement`).

![Spearman Correlation Mapping](/Users/guardian/.gemini/antigravity/brain/a7afde04-87b0-4744-b59e-ea6493409fed/artifacts/correlation_heatmap.png)

### Diabetes Class Balance
A visualization to quickly parse the volume of the dataset belonging to distinct study groups (Healthy, Type 1/2, Prediabetic).

![Diabetes Type Imbalance Summary](/Users/guardian/.gemini/antigravity/brain/a7afde04-87b0-4744-b59e-ea6493409fed/artifacts/diabetes_class_balance.png)

---

## Brainstorming: Cognition & Subjective Management Proxies

During execution, the script also outlined subjective ways to quantify cognitive deficits utilizing routine diabetes management as a proxy marker. 
> [!TIP]
> 1. **Peripheral Neuropathy Monofilament Testing:** Diminished physical sensation severely limits sustained physical activity routines (walking/exercise), concurrently correlating with physical distress and advanced cognitive decline.
> 2. **Continuous Glucose Monitor (CGM) Gap Frequency:** Lower cognition/executive functions may distinctly exhibit as poorer adherence to consistently wearing, uploading, or substituting active wearable devices.
> 3. **Correction Reactivity Windows:** A time-to-react proxy measuring how rapidly a patient acknowledges and treats severe hypoglycemia signals on a CGM, acting as a continuous executive response test.

All output datasets are located within the newly created `./eda_advanced_outputs` directory.
