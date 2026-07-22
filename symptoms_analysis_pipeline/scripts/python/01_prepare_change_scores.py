from __future__ import annotations

from pathlib import Path
import sys
import yaml
import pandas as pd


def load_config(root: Path) -> dict:
    with open(root / "config" / "config.yaml", "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    cfg = load_config(root)

    input_file = root / cfg["paths"]["input_excel"]
    sheet_name = cfg["paths"]["input_sheet"]
    out_dir = root / cfg["paths"]["processed_dir"]
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / cfg["files"]["per_patient_differences"]

    participant_col = cfg["columns"]["participant"]
    group_col = cfg["columns"]["group"]
    timepoint_col = cfg["columns"]["timepoint"]
    scores = list(cfg["columns"]["scores"])

    baseline = cfg["timepoints"]["baseline"]
    followups = [cfg["timepoints"]["end_part_1"], cfg["timepoints"]["end_part_2"], cfg["timepoints"]["end_part_3"]]

    required_cols = [participant_col, group_col, timepoint_col, *scores]

    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    df = pd.read_excel(input_file, sheet_name=sheet_name)

    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = df[required_cols].copy()
    for c in scores:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    wanted_tps = [baseline, *followups]
    df = df[df[timepoint_col].isin(wanted_tps)].copy()

    wide = (
        df.groupby([participant_col, group_col, timepoint_col], as_index=False)[scores]
        .mean()
        .pivot(index=[participant_col, group_col], columns=timepoint_col, values=scores)
    )

    pieces = []
    for followup in followups:
        if followup not in wide.columns.get_level_values(1):
            continue
        if baseline not in wide.columns.get_level_values(1):
            continue

        follow_vals = wide.xs(followup, axis=1, level=1)
        base_vals = wide.xs(baseline, axis=1, level=1)
        diff = base_vals.reindex(columns=scores) - follow_vals.reindex(columns=scores)
        diff = diff.reset_index()
        diff.insert(2, "tp_from", followup)
        diff.insert(3, "tp_to", baseline)
        pieces.append(diff)

    if not pieces:
        raise ValueError("No follow-up timepoints were found in the input data.")

    out = pd.concat(pieces, ignore_index=True)
    out = out[[participant_col, group_col, "tp_from", "tp_to", *scores]]
    out.to_excel(out_file, index=False)
    print(f"Wrote: {out_file}")
    print(out.head())


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise
