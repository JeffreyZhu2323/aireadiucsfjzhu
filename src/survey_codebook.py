"""
Codebook for the five AI-READI participant surveys held in observation.csv.

Answers in the OMOP observation table are stored as numeric codes in
value_as_number; observation_source_value holds "<item_code>, <truncated item
text>". This module is the single source of truth for item membership, item
text, and code -> label mappings, and is imported by
add_surveys_to_full_person_aggregates.py.

How the code -> label mappings were established (no concept table ships with
this export, so mappings were recovered from the survey PDFs plus the
value_as_concept_id column):

  1. AI-READI local concept ids (2005200xxx) were assigned in survey display
     order. dmlmd maps 10 -> 2005200332, 5 -> 2005200333, 0 -> 2005200334,
     which fixes 10 as the first (best) option for every dml* item.
  2. Standard concepts shared across items pin down shared answer text.
     45876662 is used by pxahc1=0 and pxdhc1..7=1, so it is "Never".
     2005200238 is used by diet1=0 and diet5=2, and the only answer text those
     two options share is "Less than 1 time", so diet code 0 is the first
     (healthiest) option. 46238128 is shared by diet1/diet6/diet7=2, all
     "4 or more times", confirming code 2 is the last option.
     45882470 is used by pxji1=0 and pxji2=0; "Don't know" is the only option
     those two questions have in common, so code 0 is a MISSING code there.
  3. Sentinels: 36311156=777 "Prefer not to say", 45876636=555 "Not
     applicable / don't know", 45884425=999 "Not sure".
  4. Ordinal direction was cross-checked against marginal distributions and,
     for pxed1, against the independent years_of_education observation item
     (codes 2..11 line up with self-reported grade completed).

Items whose mapping is inferred rather than verified are listed in
INFERRED_MAPPINGS and are deliberately excluded from derived scores.
"""

from __future__ import annotations

# Codes that mean "no substantive answer" across every survey item here.
MISSING_CODES = (555.0, 777.0, 999.0)

# Items where 0 additionally means "Don't know" rather than a scale point.
ZERO_IS_MISSING = ("pxji1", "pxji2")

PNS = "Prefer not to say"
NA_DK = "Not applicable / don't know"
NOT_SURE = "Not sure"

_SENTINELS = {555.0: NA_DK, 777.0: PNS, 999.0: NOT_SURE}

YES_NO = {1.0: "Yes", 0.0: "No"}

# --------------------------------------------------------------------------
# CES-D-10 depression (items scored 0-3; items 5 and 8 are already
# reverse-scored at capture, verified by their positive correlation with
# ces3 "I felt depressed": +0.31 and +0.49 Spearman).
# --------------------------------------------------------------------------
CESD_FREQ = {
    0.0: "Rarely or none of the time (<1 day)",
    1.0: "Some or a little of the time (1-2 days)",
    2.0: "Occasionally or a moderate amount of time (3-4 days)",
    3.0: "All of the time (5-7 days)",
}

CESD_ITEMS = {
    "ces1": "I was bothered by things that usually don't bother me",
    "ces2": "I had trouble keeping my mind on what I was doing",
    "ces3": "I felt depressed",
    "ces4": "I felt that everything I did was an effort",
    "ces5": "I felt hopeful about the future (reverse-scored)",
    "ces6": "I felt fearful",
    "ces7": "My sleep was restless",
    "ces8": "I was happy (reverse-scored)",
    "ces9": "I felt lonely",
    "ces10": "I could not 'get going'",
}

# Validated CES-D-10 three-factor structure.
CESD_DEPRESSED_AFFECT = ["ces1", "ces3", "ces6", "ces9"]
CESD_SOMATIC = ["ces2", "ces4", "ces7", "ces10"]
CESD_ANHEDONIA = ["ces5", "ces8"]
CESD_CUTOFF = 10  # score >= 10 screens positive for depression (per survey PDF)

# --------------------------------------------------------------------------
# PAID-5 diabetes distress (items scored 0-4).
# --------------------------------------------------------------------------
PAID_SEVERITY = {
    0.0: "Not a problem",
    1.0: "Minor problem",
    2.0: "Moderate problem",
    3.0: "Somewhat serious problem",
    4.0: "Serious problem",
}

PAID_ITEMS = {
    "paid_scrd": "Feeling scared when you think about living with diabetes",
    "paid_dpr": "Feeling depressed when you think about living with diabetes",
    "paid_wr": "Worrying about the future and the possibility of serious complications",
    "paid_eng": "Feeling that diabetes is taking up too much mental and physical energy every day",
    "paid_cml": "Coping with complications of diabetes",
}
PAID_CUTOFF = 8  # score >= 8 indicates diabetes-related emotional distress (per PDF)

# --------------------------------------------------------------------------
# Dietary questionnaire, 9 items scored 0/1/2 with 0 = healthiest option.
# Higher total = less healthy diet.
# --------------------------------------------------------------------------
DIET_ITEMS = {
    "diet1": "How many times a week did you eat fast food meals or snacks?",
    "diet2": "How many servings of fruit did you eat each day?",
    "diet3": "How many servings of vegetables did you eat each day?",
    "diet4": "How many regular sodas or glasses of sweet tea did you drink each day?",
    "diet5": "How many times a week did you eat beans, chicken, or fish?",
    "diet6": "How many times a week did you eat regular snack chips or crackers (not low-fat)?",
    "diet7": "How many times a week did you eat desserts and other sweets (not low-fat)?",
    "diet8": "How much margarine, butter, or meat fat do you use to season vegetables or put on potatoes, bread or corn?",
    "diet9": "How many servings of fruit juice did you drink each day?",
}

DIET_LABELS = {
    "diet1": {0.0: "Less than 1 time", 1.0: "1-3 times", 2.0: "4 or more times"},
    "diet2": {0.0: "5 or more", 1.0: "3-4", 2.0: "2 or less"},
    "diet3": {0.0: "5 or more", 1.0: "3-4", 2.0: "2 or less"},
    "diet4": {0.0: "Less than 1", 1.0: "1-2", 2.0: "3 or more"},
    "diet5": {0.0: "3 or more times", 1.0: "1-2 times", 2.0: "Less than 1 time"},
    "diet6": {0.0: "1 time or less", 1.0: "2-3 times", 2.0: "4 or more times"},
    "diet7": {0.0: "1 time or less", 1.0: "2-3 times", 2.0: "4 or more times"},
    "diet8": {0.0: "Very little", 1.0: "Some", 2.0: "A lot"},
    "diet9": {0.0: "Less than 1", 1.0: "1-2", 2.0: "3 or more"},
}

DIET_SUGARY_DRINKS = ["diet4", "diet9"]           # sodas/sweet tea + fruit juice
DIET_FRUIT_VEG = ["diet2", "diet3"]               # higher = fewer servings
DIET_FAT_PROCESSED = ["diet1", "diet6", "diet7", "diet8"]

# --------------------------------------------------------------------------
# Diabetes self-care, 10 items scored 10/5/0 with 10 = best self-care.
# Higher total = better self-care. Domains follow the survey's own framing
# ("dietary habits, exercise, healthcare, and foot care").
# --------------------------------------------------------------------------
DM_ITEMS = {
    "dmlvex": "In a typical week, how often do you engage in vigorous exercise?",
    "dmlact": "What would you consider your typical activity level to be?",
    "dmlhex": "How often would you say that you engage in daily home exercise?",
    "dmlsugar": "How often do you consume simple sugars?",
    "dmlgrain": "In a typical day, how often do you consume whole grains?",
    "dmlpor": "Which of the following best describes your approach to portion control?",
    "dmlfrveg": "In a typical day, how many servings of fruits and vegetables do you consume?",
    "dmledu": "How often would you say that you engage in diabetes health education?",
    "dmlfeet": "How often do you inspect your feet?",
    "dmlmd": "Which option best describes you with respect to physician care and medications?",
}

DM_LABELS = {
    "dmlvex": {
        10.0: "At least 30 minutes daily 3 or more days a week",
        5.0: "Less than 30 minutes daily less than 3 days a week",
        0.0: "Rarely or none",
    },
    "dmlact": {
        10.0: "Active (stairs, walks, housework)",
        5.0: "Somewhat active",
        0.0: "Not active (TV, computer, phone, use lifts)",
    },
    "dmlhex": {
        10.0: "Doing daily or almost daily",
        5.0: "Irregular; once or twice a week",
        0.0: "Rarely or none",
    },
    "dmlsugar": {
        10.0: "Rarely eat sweets",
        5.0: "Occasionally eat sweets",
        0.0: "Frequently eat sweets",
    },
    "dmlgrain": {
        10.0: "Usually eat whole grain",
        5.0: "Rarely eat whole grain",
        0.0: "Usually eat white bread or white rice",
    },
    "dmlpor": {
        10.0: "Eat small portions and avoid second servings",
        5.0: "Occasionally limit portions",
        0.0: "Rarely or none; eat at buffets and restaurants with large portions",
    },
    "dmlfrveg": {
        10.0: "4-5 servings or pieces of fruits and raw vegetables per day",
        5.0: "1-3 servings or pieces per day",
        0.0: "Rarely or none",
    },
    "dmledu": {
        10.0: "Regularly; at least monthly",
        5.0: "Occasionally; every 2-3 months",
        0.0: "Rarely or none",
    },
    "dmlfeet": {10.0: "Daily or weekly", 5.0: "Monthly", 0.0: "Rarely or none"},
    "dmlmd": {
        10.0: "3 or more visits a year to the same doctor; regular with medicines",
        5.0: "1-2 visits a year; forgetting medicines frequently",
        0.0: "Rarely or none; not taking 1 or more recommended medicines",
    },
}

DM_EXERCISE = ["dmlvex", "dmlact", "dmlhex"]
DM_DIET = ["dmlsugar", "dmlgrain", "dmlpor", "dmlfrveg"]
DM_HEALTHCARE = ["dmledu", "dmlmd"]
DM_FOOTCARE = ["dmlfeet"]

# --------------------------------------------------------------------------
# PhenX Social Determinants of Health (36 items, n=2063).
# --------------------------------------------------------------------------
_TIME_SINCE_VISIT = {
    0.0: "Never",
    1.0: "Within the past year",
    2.0: "Within the last 2 years",
    3.0: "Within the last 3 years",
    4.0: "Within the last 5 years",
    5.0: "Within the last 10 years",
    6.0: "10 years ago or more",
}

_VISIT_COUNT = {
    # Spelled out rather than "None" so the label survives a CSV round-trip,
    # where pandas would otherwise read it back as a missing value.
    0.0: "None (0 visits)",
    1.0: "1-10 times",
    2.0: "11-19 times",
    3.0: "20-29 times",
    4.0: "30-39 times",
    5.0: "40-49 times",
    6.0: "50 or more times",
}

_NEVER_ALWAYS = {
    1.0: "Never",
    2.0: "Rarely",
    3.0: "Sometimes",
    4.0: "Most of the time",
    5.0: "Always",
}

_COVERAGE = {1.0: "Covered", 0.0: "Not covered", 999.0: NOT_SURE}

EDUCATION_LABELS = {
    0.0: "Never attended / kindergarten only",
    1.0: "1st grade",
    2.0: "2nd grade",
    3.0: "3rd grade",
    4.0: "4th grade",
    5.0: "5th grade",
    6.0: "6th grade",
    7.0: "7th grade",
    8.0: "8th grade",
    9.0: "9th grade",
    10.0: "10th grade",
    11.0: "11th grade",
    12.0: "12th grade, no diploma",
    13.0: "High school graduate",
    14.0: "GED or equivalent",
    15.0: "Some college, no degree",
    16.0: "Associate degree: occupational, technical, or vocational",
    17.0: "Associate degree: academic",
    18.0: "Bachelor's degree",
    19.0: "Master's degree",
    20.0: "Professional school degree (MD, DDS, DVM, JD)",
    21.0: "Doctoral degree (PhD, EdD)",
}

# Conventional years-of-schooling equivalents for each attainment level.
EDUCATION_YEARS = {
    0.0: 0.0, 1.0: 1.0, 2.0: 2.0, 3.0: 3.0, 4.0: 4.0, 5.0: 5.0, 6.0: 6.0,
    7.0: 7.0, 8.0: 8.0, 9.0: 9.0, 10.0: 10.0, 11.0: 11.0,
    12.0: 12.0,   # 12th grade, no diploma
    13.0: 12.0,   # high school graduate
    14.0: 12.0,   # GED
    15.0: 13.0,   # some college
    16.0: 14.0,   # associate, vocational
    17.0: 14.0,   # associate, academic
    18.0: 16.0,   # bachelor's
    19.0: 18.0,   # master's
    20.0: 20.0,   # professional
    21.0: 21.0,   # doctoral
}
EDUCATION_COLLEGE_GRAD_MIN = 18.0  # bachelor's or higher

HOUSING_STATUS_LABELS = {
    1.0: "I have a steady place to live",
    2.0: "I have a place to live today, but I am worried about losing it",
    0.0: "I do not have a steady place to live",
}
HOUSING_UNSTABLE_CODES = (0.0, 2.0)  # worried about losing it, or no steady place

# pxhi2 is the one multi-select item: up to 6 rows per person, code 8 = none.
HOUSING_PROBLEM_LABELS = {
    1.0: "Pests such as bugs, ants, or mice",
    2.0: "Mold",
    3.0: "Lead paint or pipes",
    4.0: "Lack of heat",
    5.0: "Oven or stove not working",
    6.0: "Smoke detectors missing or not working",
    7.0: "Water leaks",
    8.0: "None of the above",
}
HOUSING_PROBLEM_SLUGS = {
    1.0: "pests",
    2.0: "mold",
    3.0: "lead",
    4.0: "no_heat",
    5.0: "broken_stove",
    6.0: "no_smoke_detector",
    7.0: "water_leaks",
}
HOUSING_PROBLEM_NONE = 8.0

USUAL_CARE_PLACE_LABELS = {
    1.0: "A doctor's office or health center",
    2.0: "Walk-in clinic, urgent care center, or retail clinic",
    3.0: "Emergency room",
    4.0: "A VA Medical Center or VA outpatient clinic",
    5.0: "Some other place",
    6.0: "Does not go to one place most often",
}

USUAL_CARE_PLACE_EXISTS_LABELS = {
    1.0: "Yes",
    0.0: "There is NO place",
    2.0: "There is MORE THAN ONE place",
}
NO_USUAL_CARE_CODE = 0.0

JOB_LOSS_LABELS = {
    0.0: "Don't know",
    1.0: "Very likely",
    2.0: "Fairly likely",
    3.0: "Not too likely",
    4.0: "Not at all likely",
    5.0: "Will be leaving the labor force (vol.)",
}
JOB_LOSS_LIKELY_CODES = (1.0, 2.0)      # very or fairly likely
JOB_LOSS_SCORABLE_CODES = (1.0, 2.0, 3.0, 4.0)

JOB_FIND_LABELS = {
    0.0: "Don't know",
    1.0: "Very easy",
    2.0: "Somewhat easy",
    3.0: "Not easy at all",
}

ENGLISH_PROFICIENCY_LABELS = {
    1.0: "Very well",
    2.0: "Well",
    3.0: "Not well",
    0.0: "Not at all",
}

SDOH_ITEMS = {
    "pxep1": "Is English your first language?",
    "pxep2": "Do you speak another language at home?",
    "pxep3": "How well do you speak English?",
    "pxed1": "What is the highest grade or level of school you have completed?",
    "pxji1": "Thinking about the next 12 months, how likely is it that you will lose your job or be laid off?",
    "pxji2": "How easy would it be to find a job with another employer with about the same income and benefits?",
    "pxhi1": "What is your living situation today?",
    "pxhi2": "Do you have problems with any of the following where you live? (choose all that apply)",
    "pxpa1": "In the last 12 months, have you delayed getting prescription medicines because of worry about the cost?",
    "pxpa2": "Did any other person in your household delay getting prescription medicines because of cost?",
    "pxpa3": "In the last 12 months, did you need prescription medicines but not get them because you couldn't afford it?",
    "pxpa4": "In the last 12 months, did anyone else in your household need prescription medicines but not get them because of cost?",
    "pxhic1": "Covered by insurance through a current or former employer or union?",
    "pxhic2": "Covered by insurance purchased directly from an insurance company?",
    "pxhic3": "Covered by Medicare?",
    "pxhic4": "Covered by Medicaid, Medical Assistance, CHIP, or other state/government assistance?",
    "pxhic5": "Covered by TRICARE or other military health care, including VA?",
    "pxhic7": "Covered by any other type of health insurance or health coverage plan?",
    "pxhic8": "Does this mean you currently have no health insurance or health coverage plan?",
    "pxahc1": "About how long has it been since you last saw a doctor or other health care professional?",
    "pxahc2": "Was this a wellness visit, physical, or general purpose check-up?",
    "pxahc3": "About how long has it been since you last saw a doctor for a wellness visit or check-up?",
    "pxahc4": "Is there a place that you USUALLY go to if you are sick and need health care?",
    "pxahc5": "What kind of place is it / do you go to most often?",
    "pxahc6": "During the past 12 months, how many times have you gone to a hospital emergency room?",
    "pxahc7": "During the past 12 months, how many times have you gone to an urgent care center or retail clinic?",
    "pxahc8": "During the past 12 months, have you been hospitalized overnight?",
    "pxahc9": "During the past 12 months, have you DELAYED getting medical care because of the cost?",
    "pxahc10": "During the past 12 months, did you need medical care but NOT GET IT because of the cost?",
    "pxdhc1": "You are treated with less courtesy than other people.",
    "pxdhc2": "You are treated with less respect than other people.",
    "pxdhc3": "You receive poorer service than others.",
    "pxdhc4": "A doctor or nurse acts as if he or she thinks you are not smart.",
    "pxdhc5": "A doctor or nurse acts as if he or she is afraid of you.",
    "pxdhc6": "A doctor or nurse acts as if he or she is better than you.",
    "pxdhc7": "You feel like a doctor or nurse is not listening to what you were saying.",
}

SDOH_LABELS = {
    "pxep1": YES_NO,
    "pxep2": YES_NO,
    "pxep3": ENGLISH_PROFICIENCY_LABELS,
    "pxed1": EDUCATION_LABELS,
    "pxji1": JOB_LOSS_LABELS,
    "pxji2": JOB_FIND_LABELS,
    "pxhi1": HOUSING_STATUS_LABELS,
    "pxhi2": HOUSING_PROBLEM_LABELS,
    "pxpa1": YES_NO,
    "pxpa2": YES_NO,
    "pxpa3": YES_NO,
    "pxpa4": YES_NO,
    "pxhic1": _COVERAGE,
    "pxhic2": _COVERAGE,
    "pxhic3": _COVERAGE,
    "pxhic4": _COVERAGE,
    "pxhic5": _COVERAGE,
    "pxhic7": _COVERAGE,
    "pxhic8": {},  # ambiguous, see INFERRED_MAPPINGS
    "pxahc1": _TIME_SINCE_VISIT,
    "pxahc2": YES_NO,
    "pxahc3": _TIME_SINCE_VISIT,
    "pxahc4": USUAL_CARE_PLACE_EXISTS_LABELS,
    "pxahc5": USUAL_CARE_PLACE_LABELS,
    "pxahc6": _VISIT_COUNT,
    "pxahc7": _VISIT_COUNT,
    "pxahc8": YES_NO,
    "pxahc9": YES_NO,
    "pxahc10": YES_NO,
    "pxdhc1": _NEVER_ALWAYS,
    "pxdhc2": _NEVER_ALWAYS,
    "pxdhc3": _NEVER_ALWAYS,
    "pxdhc4": _NEVER_ALWAYS,
    "pxdhc5": _NEVER_ALWAYS,
    "pxdhc6": _NEVER_ALWAYS,
    "pxdhc7": _NEVER_ALWAYS,
}

INSURANCE_ITEMS = ["pxhic1", "pxhic2", "pxhic3", "pxhic4", "pxhic5", "pxhic7"]
RX_COST_ITEMS = ["pxpa1", "pxpa2", "pxpa3", "pxpa4"]
RX_COST_SELF_ITEMS = ["pxpa1", "pxpa3"]
MED_DISCRIMINATION_ITEMS = [f"pxdhc{i}" for i in range(1, 8)]

# --------------------------------------------------------------------------
# PhenX Food Insecurity (USDA-style, n=2280).
# pxfi1/pxfi2 are the two Hunger Vital Sign screener items. Their code 0 is
# "Never true"; codes 1 and 2 are "Often true" and "Sometimes true" but which
# is which could not be verified, so both are treated as affirmative -- the
# standard scoring rule anyway, which makes the ambiguity moot.
# --------------------------------------------------------------------------
FOOD_INSECURITY_ITEMS = {
    "pxfi1": "The food that (I/we) bought just didn't last, and (I/we) didn't have money to get more",
    "pxfi2": "(I/we) couldn't afford to eat balanced meals",
    "pxfi3": "In the last year, did you and/or other adults cut the size of meals or skip meals because there wasn't enough money for food?",
    "pxfi4": "In the last 12 months, did you ever eat less than you felt you should because there wasn't enough money for food?",
    "pxfi5": "In the last 12 months, were you ever hungry but didn't eat because there wasn't enough money for food?",
    "food_insecurity_change_diet_frequency": "How often did this happen? (follow-up to pxfi3)",
}

FI_OFTEN_SOMETIMES = {0.0: "Never true", 1.0: "Often or sometimes true (a)", 2.0: "Often or sometimes true (b)"}

FOOD_INSECURITY_LABELS = {
    "pxfi1": FI_OFTEN_SOMETIMES,
    "pxfi2": FI_OFTEN_SOMETIMES,
    "pxfi3": YES_NO,
    "pxfi4": YES_NO,
    "pxfi5": YES_NO,
    "food_insecurity_change_diet_frequency": {},  # direction unverified
}

# Affirmative (food-insecure) response codes per item.
FI_AFFIRMATIVE = {
    "pxfi1": (1.0, 2.0),
    "pxfi2": (1.0, 2.0),
    "pxfi3": (1.0,),
    "pxfi4": (1.0,),
    "pxfi5": (1.0,),
}
FI_SCORE_ITEMS = ["pxfi1", "pxfi2", "pxfi3", "pxfi4", "pxfi5"]
FI_HVS_ITEMS = ["pxfi1", "pxfi2"]
FI_CUTOFF = 2  # >= 2 affirmatives indicates food insecurity (USDA short form)

# --------------------------------------------------------------------------
# Mappings that could not be fully verified. Excluded from derived scores;
# raw codes are still written to the item-level file.
# --------------------------------------------------------------------------
INFERRED_MAPPINGS = {
    "pxep3": (
        "Code 1/2/3 = Very well/Well/Not well fits the marginal counts "
        "(174/143/12) and the contiguous concept ids 45877707/05/06; code 0 "
        "(n=2) is assigned 'Not at all' by elimination. Used for labels only, "
        "never for a derived score."
    ),
    "pxfi1/pxfi2": (
        "Codes 1 and 2 are 'Often true' and 'Sometimes true' in unknown order. "
        "Both count as affirmative under standard scoring, so scores are "
        "unaffected; the labels say 'Often or sometimes true (a)/(b)'."
    ),
    "food_insecurity_change_diet_frequency": (
        "USDA frequency follow-up (almost every month / some months / only 1-2 "
        "months) in unknown code order, so it is EXCLUDED from "
        "food_insecurity_score. That makes the score a 5-item affirmative "
        "count rather than the USDA 6-item raw score."
    ),
    "pxhic8": (
        "Only 31 respondents and the two codes cannot be resolved against the "
        "'I do NOT have insurance' / 'I HAVE some kind of insurance' options. "
        "Not used; sdoh_uninsured is derived from pxhic1..7 instead."
    ),
}

ALL_ITEM_TEXT = {
    **CESD_ITEMS,
    **PAID_ITEMS,
    **DIET_ITEMS,
    **DM_ITEMS,
    **SDOH_ITEMS,
    **FOOD_INSECURITY_ITEMS,
}

ALL_LABELS = {
    **{k: CESD_FREQ for k in CESD_ITEMS},
    **{k: PAID_SEVERITY for k in PAID_ITEMS},
    **DIET_LABELS,
    **DM_LABELS,
    **SDOH_LABELS,
    **FOOD_INSECURITY_LABELS,
}

SURVEY_OF_ITEM = {
    **{k: "CES-D-10 depression" for k in CESD_ITEMS},
    **{k: "PAID-5 diabetes distress" for k in PAID_ITEMS},
    **{k: "Dietary questionnaire" for k in DIET_ITEMS},
    **{k: "Diabetes self-care" for k in DM_ITEMS},
    **{k: "PhenX SDOH" for k in SDOH_ITEMS},
    **{k: "PhenX food insecurity" for k in FOOD_INSECURITY_ITEMS},
}

# Nominal (unordered) items that get a companion *_label column in the
# item-level export; everything else is ordinal and readable as a number.
NOMINAL_ITEMS = [
    "pxed1", "pxhi1", "pxji1", "pxji2", "pxep3",
    "pxahc1", "pxahc3", "pxahc4", "pxahc5", "pxahc6", "pxahc7",
    "pxfi1", "pxfi2",
]

ALL_ITEMS = list(ALL_ITEM_TEXT)


def sentinel_label(code: float) -> str | None:
    """Human-readable name for a missing-data sentinel code."""
    return _SENTINELS.get(code)
