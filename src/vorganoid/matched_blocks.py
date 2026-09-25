"""Donor-equal, plate-and-forskolin-dose-matched FIS size sensitivity.

The included OrgaSegment thresholds were chosen after looking at that same accession.
This computes an assay sensitivity result, not a treatment recommendation or diagnosis.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.stats import binomtest


def matched_block_effects(df: pd.DataFrame, drug: str, control: str = "DMSO",
                          small_max: float = 722., large_min: float = 1400.,
                          min_per_cell: int = 3) -> dict:
    required = {"donor", "experiment", "forskolin_concentration_µM", "condition", "A0", "swelling"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"missing columns: {missing}")
    if not (0 < small_max < large_min and min_per_cell >= 1):
        raise ValueError("require 0 < small_max < large_min and min_per_cell >= 1")
    y = df[df.condition.isin([drug, control])].copy()
    y = y[(y["forskolin_concentration_µM"] > 0) & (y.A0 > 0) & (y.swelling > 0)]
    y["size"] = np.select([y.A0 < small_max, y.A0 >= large_min], ["small", "large"], default="middle")
    y = y[y["size"] != "middle"].copy()
    y["log_swelling"] = np.log(y.swelling)
    cells = [(c, size) for c in (control, drug) for size in ("small", "large")]
    blocks = []
    for (donor, experiment, dose), z in y.groupby(["donor", "experiment", "forskolin_concentration_µM"]):
        q = z.groupby(["condition", "size"]).log_swelling.agg(["mean", "size"])
        if any(cell not in q.index or q.loc[cell, "size"] < min_per_cell for cell in cells):
            continue
        effect = (q.loc[(drug, "large"), "mean"] - q.loc[(control, "large"), "mean"]
                  - q.loc[(drug, "small"), "mean"] + q.loc[(control, "small"), "mean"])
        blocks.append({"donor": str(donor), "experiment": str(experiment), "dose_uM": float(dose),
                       "effect": float(effect),
                       "cell_counts": {f"{c}:{size}": int(q.loc[(c, size), "size"]) for c, size in cells}})
    donors = []
    for donor, z in pd.DataFrame(blocks).groupby("donor") if blocks else []:
        donors.append({"donor": donor, "n_blocks": len(z), "median_effect": float(z.effect.median())})
    n, positive = len(donors), sum(d["median_effect"] > 0 for d in donors)
    p = float(binomtest(positive, n, .5, alternative="greater").pvalue) if n else None
    return {"drug": drug, "control": control, "small_max": small_max, "large_min": large_min,
            "min_per_cell": min_per_cell, "n_donors": n, "n_positive_donors": positive,
            "one_sided_sign_p": p, "verdict": "UNINTERPRETABLE" if n < 8 else
            "PASS" if positive > n / 2 and p < .05 else "FAIL",
            "warning": "same-source post-hoc size thresholds; donor sign test is not independent replication or clinical validation",
            "donors": donors, "blocks": blocks}
