from pathlib import Path
import pandas as pd

CONTAINER = "jeffreyzhuaireadi"
UUID = "7bc8c89b-8457-41fc-917c-b5571ff90ea0"

ROOT = Path(__file__).resolve().parents[1]
MAN = ROOT / "data" / "manifests" / CONTAINER / UUID / "dataset"
OUTDIR = ROOT / "outputs"
LOGDIR = ROOT / "logs"
BASE_DIR = Path(__file__).resolve().parent.parent
MANIFEST_DIR = BASE_DIR / "data/manifests/jeffreyzhuaireadi/7bc8c89b-8457-41fc-917c-b5571ff90ea0/dataset/"
RAW_BASE = BASE_DIR / "data/raw/jeffreyzhuaireadi/jeffreyzhuaireadi/7bc8c89b-8457-41fc-917c-b5571ff90ea0/dataset/"
OUT_GA = OUTDIR / "manifest_lists" / "files_glucose_activity.txt"
OUT_ENV = OUTDIR / "manifest_lists" / "files_environment.txt"

OUTDIR.mkdir(parents=True, exist_ok=True)
LOGDIR.mkdir(parents=True, exist_ok=True)
OUT_GA.parent.mkdir(parents=True, exist_ok=True)


def _norm(p: str) -> str:
    """
    Normalize manifest file paths to include-path format relative to the container root.
    """
    p = str(p).lstrip("/")
    if not p.startswith(f"{UUID}/dataset/"):
        p = f"{UUID}/dataset/{p}"
    return p


def write_manifest_file_lists() -> None:
    """
    Read manifests and write out text files listing all raw glucose/activity and env paths.

    This is intended as a one-off utility and is NOT run on import; call explicitly or
    run this module as a script.
    """
    glucose = pd.read_csv(MAN / "wearable_blood_glucose" / "manifest.tsv", sep="\t")
    activity = pd.read_csv(MAN / "wearable_activity_monitor" / "manifest.tsv", sep="\t")
    env = pd.read_csv(MAN / "environment" / "manifest.tsv", sep="\t")

    paths_ga: set[str] = set()
    paths_env: set[str] = set()

    # glucose raw files
    if "glucose_filepath" in glucose.columns:
        for p in glucose["glucose_filepath"].dropna().astype(str):
            paths_ga.add(_norm(p))

    # activity raw files (ALL *_filepath columns)
    filepath_cols = [c for c in activity.columns if c.endswith("_filepath")]
    for c in filepath_cols:
        for p in activity[c].dropna().astype(str):
            paths_ga.add(_norm(p))

    # environment raw files
    if "env_sensor_filepath" in env.columns:
        for p in env["env_sensor_filepath"].dropna().astype(str):
            paths_env.add(_norm(p))

    OUT_GA.write_text("\n".join(sorted(paths_ga)) + "\n")
    OUT_ENV.write_text("\n".join(sorted(paths_env)) + "\n")

    print("wrote:", OUT_GA, "n=", len(paths_ga))
    print("wrote:", OUT_ENV, "n=", len(paths_env))


def main() -> None:
    """CLI entrypoint for generating manifest-based file lists."""
    write_manifest_file_lists()


if __name__ == "__main__":
    main()