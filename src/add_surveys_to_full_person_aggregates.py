"""
Extract numeric and categorical features from the five participant surveys in
observation.csv and add them to outputs/full_person_aggregates.csv.

Surveys covered (item membership and code -> label maps live in
src/survey_codebook.py):
  * CES-D-10 depression            ces1..ces10
  * PAID-5 diabetes distress       paid_scrd, paid_dpr, paid_wr, paid_eng, paid_cml
  * Dietary questionnaire          diet1..diet9
  * Diabetes self-care             10 dml* items
  * PhenX SDOH                     36 items, 8 domains
  * PhenX food insecurity          pxfi1..pxfi5 (+ frequency follow-up)

Two scoring decisions worth knowing about:

  * dietscore in the export is unusable. For 1435 of 2151 respondents it is
    stored as exactly (true score - 20) as a negative number in
    value_as_string, with value_as_number null; the other 716 are correct.
    diet_risk_score here is recomputed from diet1..diet9 and reproduces the
    true 0-18 score for everyone.
  * cestl and paidscore ARE exact sums of their items (verified for all
    2224 and 2207 complete respondents), but they are recomputed anyway so
    that the proration rule below applies uniformly.

Missing data: codes 555/777/999 (not applicable, prefer not to say, not sure)
become NaN, as does code 0 on pxji1/pxji2 where it means "Don't know". A scale
score is emitted when at least 80% of its items were answered, prorated to the
full item count; otherwise NaN. For the 10-item CES-D-10 that rule permits at
most 2 missing items, matching the instrument's own scoring instruction.

Writes three files:
  * outputs/full_person_aggregates.csv     updated in place, survey block
                                           inserted after moca_score
  * outputs/survey_items_by_person.csv     one row per person, one column per
                                           survey item (cleaned codes)
  * outputs/survey_codebook.csv            tidy code -> label codebook with
                                           observed response counts

Run from repo root:
    .venv/bin/python src/add_surveys_to_full_person_aggregates.py
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

import survey_codebook as cb

ROOT = Path(__file__).resolve().parents[1]
OBSERVATION = ROOT / "data" / "raw" / "omop" / "observation.csv"
FULL_PATH = ROOT / "outputs" / "cohort" / "full_person_aggregates.csv"
ITEMS_PATH = ROOT / "outputs" / "cohort" / "survey_items_by_person.csv"
CODEBOOK_PATH = ROOT / "outputs" / "cohort" / "survey_codebook.csv"

AFTER_COL = "moca_score"  # survey block is inserted after this column
MIN_ANSWERED_FRACTION = 0.8


# ---------------------------------------------------------------- extraction

def load_items() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (wide single-response items, raw rows for the multi-select item).

    The wide frame is person_id x item_code holding cleaned numeric codes, with
    every missing-data sentinel already converted to NaN.
    """
    obs = pd.read_csv(
        OBSERVATION,
        usecols=["person_id", "observation_source_value", "value_as_number"],
        low_memory=False,
    )
    obs["code"] = obs["observation_source_value"].str.split(",").str[0].str.strip()
    obs = obs[obs["code"].isin(cb.ALL_ITEMS)].copy()
    obs["person_id"] = pd.to_numeric(obs["person_id"], errors="coerce").astype("Int64")
    obs["value_as_number"] = pd.to_numeric(obs["value_as_number"], errors="coerce")

    multi = obs[obs["code"] == "pxhi2"].copy()
    single = obs[obs["code"] != "pxhi2"]

    dupes = single.groupby(["code", "person_id"]).size()
    assert dupes.max() == 1, f"unexpected repeated responses: {dupes[dupes > 1]}"

    wide = single.pivot(index="person_id", columns="code", values="value_as_number")
    wide = wide.replace(dict.fromkeys(cb.MISSING_CODES, np.nan))
    for item in cb.ZERO_IS_MISSING:
        if item in wide.columns:
            wide[item] = wide[item].replace(0.0, np.nan)

    return wide, multi


def scale_score(
    wide: pd.DataFrame, items: list[str], scale_max: float
) -> tuple[pd.Series, pd.Series]:
    """Prorated sum over `items`, plus the per-person answered-item count.

    Emitted only when at least MIN_ANSWERED_FRACTION of items were answered;
    the mean of answered items is scaled up to the full item count so that
    partial responders stay on the instrument's native range.
    """
    sub = wide.reindex(columns=items)
    answered = sub.notna().sum(axis=1)
    required = math.ceil(MIN_ANSWERED_FRACTION * len(items))
    total = (sub.mean(axis=1) * len(items)).round(2)
    total = total.where(answered >= required)
    assert total.max(skipna=True) <= scale_max + 1e-9, (
        f"{items} exceeded max {scale_max}: {total.max()}"
    )
    return total, answered


def affirmative_count(
    wide: pd.DataFrame, spec: dict[str, tuple[float, ...]]
) -> tuple[pd.Series, pd.Series]:
    """Count of items answered with an affirmative code, and answered count."""
    items = list(spec)
    sub = wide.reindex(columns=items)
    answered = sub.notna().sum(axis=1)
    flags = pd.DataFrame(
        {
            item: sub[item].isin(spec[item]).astype("float").where(sub[item].notna())
            for item in items
        }
    )
    required = math.ceil(MIN_ANSWERED_FRACTION * len(items))
    count = flags.sum(axis=1).where(answered >= required)
    return count, answered


def housing_problems(multi: pd.DataFrame, index: pd.Index) -> pd.DataFrame:
    """Per-person housing problem count and one flag per problem type."""
    m = multi[~multi["value_as_number"].isin(cb.MISSING_CODES)]
    problems = m[m["value_as_number"] != cb.HOUSING_PROBLEM_NONE]

    responded = pd.Index(m["person_id"].unique())
    count = (
        problems.groupby("person_id").size().reindex(responded, fill_value=0)
    )

    out = pd.DataFrame(index=index)
    out["sdoh_housing_problem_count"] = count.reindex(index)
    for code, slug in cb.HOUSING_PROBLEM_SLUGS.items():
        has = problems.loc[problems["value_as_number"] == code, "person_id"].unique()
        flag = pd.Series(0, index=responded, dtype="float")
        flag.loc[flag.index.isin(has)] = 1
        out[f"sdoh_housing_problem_{slug}"] = flag.reindex(index)
    return out


# ------------------------------------------------------------------ features

def build_features(wide: pd.DataFrame, multi: pd.DataFrame) -> pd.DataFrame:
    f = pd.DataFrame(index=wide.index)

    # --- CES-D-10 depression -------------------------------------------
    f["cesd10_total"], f["cesd10_items_answered"] = scale_score(
        wide, list(cb.CESD_ITEMS), 30
    )
    f["cesd10_depressed"] = (f["cesd10_total"] >= cb.CESD_CUTOFF).astype("float").where(
        f["cesd10_total"].notna()
    )
    f["cesd10_depressed_affect"], _ = scale_score(wide, cb.CESD_DEPRESSED_AFFECT, 12)
    f["cesd10_somatic"], _ = scale_score(wide, cb.CESD_SOMATIC, 12)
    f["cesd10_anhedonia"], _ = scale_score(wide, cb.CESD_ANHEDONIA, 6)

    # --- PAID-5 diabetes distress --------------------------------------
    f["paid5_total"], f["paid5_items_answered"] = scale_score(
        wide, list(cb.PAID_ITEMS), 20
    )
    f["paid5_distress"] = (f["paid5_total"] >= cb.PAID_CUTOFF).astype("float").where(
        f["paid5_total"].notna()
    )

    # --- Dietary questionnaire (higher = less healthy) -----------------
    f["diet_risk_score"], f["diet_items_answered"] = scale_score(
        wide, list(cb.DIET_ITEMS), 18
    )
    f["diet_sugary_drinks"], _ = scale_score(wide, cb.DIET_SUGARY_DRINKS, 4)
    f["diet_fruit_veg_deficit"], _ = scale_score(wide, cb.DIET_FRUIT_VEG, 4)
    f["diet_fat_processed"], _ = scale_score(wide, cb.DIET_FAT_PROCESSED, 8)

    # --- Diabetes self-care (higher = better) --------------------------
    f["dm_selfcare_total"], f["dm_selfcare_items_answered"] = scale_score(
        wide, list(cb.DM_ITEMS), 100
    )
    f["dm_selfcare_exercise"], _ = scale_score(wide, cb.DM_EXERCISE, 30)
    f["dm_selfcare_diet"], _ = scale_score(wide, cb.DM_DIET, 40)
    f["dm_selfcare_healthcare"], _ = scale_score(wide, cb.DM_HEALTHCARE, 20)
    f["dm_selfcare_footcare"], _ = scale_score(wide, cb.DM_FOOTCARE, 10)

    # --- SDOH: education ------------------------------------------------
    edu = wide["pxed1"]
    f["sdoh_education_code"] = edu
    f["sdoh_education_label"] = edu.map(cb.EDUCATION_LABELS)
    f["sdoh_education_years"] = edu.map(cb.EDUCATION_YEARS)
    f["sdoh_college_graduate"] = (
        (edu >= cb.EDUCATION_COLLEGE_GRAD_MIN).astype("float").where(edu.notna())
    )

    # --- SDOH: language -------------------------------------------------
    f["sdoh_english_first_language"] = wide["pxep1"]

    # --- SDOH: housing --------------------------------------------------
    hi1 = wide["pxhi1"]
    f["sdoh_housing_status_code"] = hi1
    f["sdoh_housing_status_label"] = hi1.map(cb.HOUSING_STATUS_LABELS)
    f["sdoh_housing_instability"] = (
        hi1.isin(cb.HOUSING_UNSTABLE_CODES).astype("float").where(hi1.notna())
    )
    hp = housing_problems(multi, wide.index)
    f["sdoh_housing_problem_count"] = hp["sdoh_housing_problem_count"]
    f["sdoh_any_housing_problem"] = (
        (hp["sdoh_housing_problem_count"] > 0)
        .astype("float")
        .where(hp["sdoh_housing_problem_count"].notna())
    )

    # --- SDOH: insurance ------------------------------------------------
    ins = wide.reindex(columns=cb.INSURANCE_ITEMS)
    ins_answered = ins.notna().sum(axis=1)
    f["sdoh_insurance_type_count"] = ins.sum(axis=1).where(ins_answered >= 5)
    f["sdoh_uninsured"] = (
        (f["sdoh_insurance_type_count"] == 0)
        .astype("float")
        .where(f["sdoh_insurance_type_count"].notna())
    )

    # --- SDOH: prescription cost barriers -------------------------------
    rx_count, rx_answered = affirmative_count(
        wide, {i: (1.0,) for i in cb.RX_COST_ITEMS}
    )
    f["sdoh_rx_cost_barrier_count"] = rx_count
    rx_self = wide.reindex(columns=cb.RX_COST_SELF_ITEMS)
    f["sdoh_rx_cost_barrier_self"] = (
        rx_self.eq(1.0).any(axis=1).astype("float").where(rx_self.notna().any(axis=1))
    )

    # --- SDOH: healthcare access barriers and utilisation ---------------
    def yes_flag(item: str) -> pd.Series:
        col = wide[item]
        return col.eq(1.0).astype("float").where(col.notna())

    no_place = (
        wide["pxahc4"].eq(cb.NO_USUAL_CARE_CODE).astype("float")
        .where(wide["pxahc4"].notna())
    )
    barriers = pd.DataFrame(
        {
            "delayed_cost": yes_flag("pxahc9"),
            "unmet_cost": yes_flag("pxahc10"),
            "no_usual_place": no_place,
            "uninsured": f["sdoh_uninsured"],
        }
    )
    b_answered = barriers.notna().sum(axis=1)
    f["sdoh_care_access_barrier_count"] = barriers.sum(axis=1).where(b_answered >= 3)
    f["sdoh_no_usual_care_place"] = no_place
    f["sdoh_hospitalized_overnight"] = wide["pxahc8"]
    f["sdoh_urgent_care_visits_code"] = wide["pxahc7"]

    # --- SDOH: discrimination in medical settings ------------------------
    dhc = wide.reindex(columns=cb.MED_DISCRIMINATION_ITEMS)
    dhc_answered = dhc.notna().sum(axis=1)
    required = math.ceil(MIN_ANSWERED_FRACTION * len(cb.MED_DISCRIMINATION_ITEMS))
    f["sdoh_med_discrimination_mean"] = (
        dhc.mean(axis=1).round(3).where(dhc_answered >= required)
    )
    f["sdoh_med_discrimination_any"] = (
        (dhc > 1.0).any(axis=1).astype("float").where(dhc_answered >= required)
    )

    # --- SDOH: job insecurity (low coverage: code 0 = "Don't know") ------
    ji1 = wide["pxji1"]
    scorable = ji1.isin(cb.JOB_LOSS_SCORABLE_CODES)
    f["sdoh_job_loss_likely"] = (
        ji1.isin(cb.JOB_LOSS_LIKELY_CODES).astype("float").where(scorable)
    )

    # --- Food insecurity -------------------------------------------------
    fi_count, fi_answered = affirmative_count(wide, cb.FI_AFFIRMATIVE)
    f["food_insecurity_score"] = fi_count
    f["food_insecurity_items_answered"] = fi_answered
    f["food_insecure"] = (
        (fi_count >= cb.FI_CUTOFF).astype("float").where(fi_count.notna())
    )
    hvs = wide.reindex(columns=cb.FI_HVS_ITEMS)
    f["food_insecurity_hvs_positive"] = (
        (hvs > 0).any(axis=1).astype("float").where(hvs.notna().all(axis=1))
    )

    # --- Cumulative SDOH adversity ---------------------------------------
    # One flag per domain so that housing, which contributes two questions,
    # cannot count twice.
    adversity = pd.DataFrame(
        {
            "food": f["food_insecure"],
            "housing": pd.concat(
                [f["sdoh_housing_instability"], f["sdoh_any_housing_problem"]], axis=1
            ).max(axis=1),
            "rx_cost": f["sdoh_rx_cost_barrier_self"],
            "care_cost": barriers[["delayed_cost", "unmet_cost"]].max(axis=1),
            "access": barriers[["no_usual_place", "uninsured"]].max(axis=1),
        }
    )
    a_answered = adversity.notna().sum(axis=1)
    f["sdoh_adversity_index"] = adversity.sum(axis=1).where(a_answered >= 4)

    return f


def build_item_export(wide: pd.DataFrame, multi: pd.DataFrame) -> pd.DataFrame:
    """Person x item table of cleaned codes, with labels for nominal items."""
    ordered: list[str] = []
    for group in (
        cb.CESD_ITEMS, cb.PAID_ITEMS, cb.DIET_ITEMS, cb.DM_ITEMS,
        cb.SDOH_ITEMS, cb.FOOD_INSECURITY_ITEMS,
    ):
        ordered += [i for i in group if i in wide.columns]

    out = wide.reindex(columns=ordered).copy()
    for item in cb.NOMINAL_ITEMS:
        if item in out.columns:
            labels = cb.ALL_LABELS.get(item, {})
            out.insert(
                out.columns.get_loc(item) + 1,
                f"{item}_label",
                out[item].map(labels),
            )

    hp = housing_problems(multi, wide.index)
    out = out.join(hp)
    out.index.name = "person_id"
    return out.reset_index()


def build_codebook(wide: pd.DataFrame, multi: pd.DataFrame) -> pd.DataFrame:
    """Tidy codebook: one row per (item, response code) with observed counts."""
    raw = pd.read_csv(
        OBSERVATION,
        usecols=["observation_source_value", "value_as_number"],
        low_memory=False,
    )
    raw["code"] = raw["observation_source_value"].str.split(",").str[0].str.strip()
    raw = raw[raw["code"].isin(cb.ALL_ITEMS)]
    raw["value_as_number"] = pd.to_numeric(raw["value_as_number"], errors="coerce")

    rows = []
    for item, grp in raw.groupby("code", sort=False):
        counts = grp["value_as_number"].value_counts().sort_index()
        labels = cb.ALL_LABELS.get(item, {})
        for value, n in counts.items():
            sentinel = cb.sentinel_label(value)
            if sentinel is not None:
                label, kind = sentinel, "missing"
            elif item in cb.ZERO_IS_MISSING and value == 0.0:
                label, kind = "Don't know", "missing"
            else:
                label, kind = labels.get(value, ""), "response"
            rows.append(
                {
                    "survey": cb.SURVEY_OF_ITEM[item],
                    "item_code": item,
                    "item_text": cb.ALL_ITEM_TEXT[item],
                    "value_code": value,
                    "value_label": label,
                    "value_kind": kind,
                    "n_respondents": int(n),
                }
            )

    order = [
        "CES-D-10 depression", "PAID-5 diabetes distress", "Dietary questionnaire",
        "Diabetes self-care", "PhenX SDOH", "PhenX food insecurity",
    ]
    book = pd.DataFrame(rows)
    book["survey"] = pd.Categorical(book["survey"], categories=order, ordered=True)
    return book.sort_values(["survey", "item_code", "value_code"], kind="stable")


# ---------------------------------------------------------------------- main

def main() -> None:
    print("Loading survey items from observation table...")
    wide, multi = load_items()
    print(f"  items: {wide.shape[1] + 1}, respondents: {wide.shape[0]}")

    features = build_features(wide, multi)
    print(f"  built {features.shape[1]} survey features")

    items = build_item_export(wide, multi)
    items.to_csv(ITEMS_PATH, index=False)
    print(f"Wrote: {ITEMS_PATH}  (rows={len(items)}, cols={len(items.columns)})")

    book = build_codebook(wide, multi)
    book.to_csv(CODEBOOK_PATH, index=False)
    print(f"Wrote: {CODEBOOK_PATH}  (rows={len(book)})")

    full = pd.read_csv(FULL_PATH)
    full["person_id"] = pd.to_numeric(full["person_id"], errors="coerce").astype("Int64")
    n0 = len(full)

    new_cols = list(features.columns)
    keep = [c for c in full.columns if c not in new_cols]
    full = full[keep]

    merged = full.merge(
        features.reset_index().rename(columns={"index": "person_id"}),
        on="person_id",
        how="left",
    )
    assert len(merged) == n0, f"row count changed: {n0} -> {len(merged)}"

    ordered: list[str] = []
    for c in keep:
        ordered.append(c)
        if c == AFTER_COL:
            ordered += new_cols
    ordered += [c for c in new_cols if c not in ordered]
    merged = merged[ordered]
    merged.to_csv(FULL_PATH, index=False)

    print(f"Wrote: {FULL_PATH}  (rows={len(merged)}, cols={len(merged.columns)})")
    print("\nCoverage of new survey features:")
    cov = pd.DataFrame(
        {
            "non_null": merged[new_cols].notna().sum(),
            "pct": (merged[new_cols].notna().mean() * 100).round(1),
        }
    )
    numeric = [c for c in new_cols if pd.api.types.is_numeric_dtype(merged[c])]
    cov["min"] = merged[numeric].min()
    cov["max"] = merged[numeric].max()
    cov["mean"] = merged[numeric].mean().round(2)
    print(cov.to_string())


if __name__ == "__main__":
    main()
