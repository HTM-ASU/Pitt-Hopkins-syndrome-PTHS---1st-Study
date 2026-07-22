from __future__ import annotations

from pathlib import Path
import sys
import yaml
import numpy as np
import pandas as pd
from scipy import stats


def load_config(root: Path) -> dict:
    with open(root / "config" / "config.yaml", "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def cohens_d_independent(x: pd.Series, y: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce").dropna()
    y = pd.to_numeric(y, errors="coerce").dropna()
    n1, n2 = len(x), len(y)
    if n1 < 2 or n2 < 2:
        return np.nan
    s1 = x.std(ddof=1)
    s2 = y.std(ddof=1)
    pooled = np.sqrt(((n1 - 1) * s1**2 + (n2 - 1) * s2**2) / (n1 + n2 - 2))
    if pooled == 0 or np.isnan(pooled):
        return np.nan
    return (x.mean() - y.mean()) / pooled


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    cfg = load_config(root)

    in_file = root / cfg["paths"]["processed_dir"] / cfg["files"]["per_patient_differences"]
    out_file = root / cfg["paths"]["processed_dir"] / cfg["files"]["ttest_results"]

    participant_col = cfg["columns"]["participant"]
    group_col = cfg["columns"]["group"]
    scores = list(cfg["columns"]["scores"])
    treated = cfg["groups"]["treated"]
    control = cfg["groups"]["control"]
    followups = [cfg["timepoints"]["end_part_1"], cfg["timepoints"]["end_part_2"], cfg["timepoints"]["end_part_3"]]

    if not in_file.exists():
        raise FileNotFoundError(f"Missing input file: {in_file}")

    df = pd.read_excel(in_file)

    required = [participant_col, group_col, "tp_from", "tp_to", *scores]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in change-score file: {missing}")

    results = []
    for tp in followups:
        for feature in scores:
            g_treated = df.loc[(df[group_col] == treated) & (df["tp_from"] == tp), feature].dropna()
            g_control = df.loc[(df[group_col] == control) & (df["tp_from"] == tp), feature].dropna()

            n1, n2 = len(g_treated), len(g_control)
            if n1 == 0 or n2 == 0:
                results.append([feature, tp, np.nan, np.nan, np.nan, np.nan, n1, n2])
                continue

            t_stat, p_value = stats.ttest_ind(g_treated, g_control, equal_var=False, alternative="greater")
            d = cohens_d_independent(g_treated, g_control)
            results.append([feature, tp, p_value, t_stat, d, abs(d) if pd.notna(d) else np.nan, n1, n2])

    results_df = pd.DataFrame(
        results,
        columns=["feature", "timepoint", "p_value", "t_statistic", "cohens_d", "abs_cohens_d", "n_treated", "n_control"],
    )

    results_df.to_excel(out_file, index=False)
    print(f"Wrote: {out_file}")
    print(results_df)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise
