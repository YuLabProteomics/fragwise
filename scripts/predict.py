#!/usr/bin/env python3
"""
Predict the optimal MS2 fragmentation setting for a list of peptides.

Usage:
  python scripts/predict.py \
      --input example/example_input.csv \
      --output example/example_output.csv

Input CSV:
  sequence    Peptide sequence (oxidized Met: M*, M(ox), or M[+16])
  charge      Precursor charge state
  peptide_id  Optional identifier

Requirements:
  pip install -r requirements.txt
"""

import argparse
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
from pyteomics import mass as pmass

# Constants 
PROTON = 1.007276466812
TMTPRO_MASS = 295.189591      
OXM_MASS = 15.994915           # oxidation mass
CYS_MASS = 57.0215             # Cys carbamidomethylation 

ALL_SETTINGS = [
    "CID_22", "CID_25", "CID_28", "CID_31", "CID_34",
    "HCD_22", "HCD_25", "HCD_28", "HCD_31", "HCD_34",
]

MODEL_TAG = "fragwise"
N_FOLDS = 5

FEATURE_ORDER = [
    "charge", "precursor_mz", "frac_H", "oxM_count",
    "oxM_position_mean", "RK_frac", "pI", "GRAVY",
    "proline_position_mean",
]

# Sequence parsing
OXM_MASS_TOL = 0.05  # Da

def _sub_oxM_bracket(match):
    """Convert M[+mass] to 'm' when mass matches oxidized Met."""
    delta = float(match.group(1))
    if abs(delta - OXM_MASS) <= OXM_MASS_TOL:
        return "m"
    return match.group(0)

def parse_sequence(raw_seq):
    """Parse a peptide sequence and identify oxidized methionines.
    Returns the clean sequence, oxidized Met count, and positions.
    Supports M*, M(ox), and M[+16] notation.
    """
    s = raw_seq.strip()

    s = re.sub(r"M\*", "m", s)
    s = re.sub(r"M\(ox\)", "m", s, flags=re.IGNORECASE)
    s = re.sub(r"M\[\+([\d.]+)\]", _sub_oxM_bracket, s)

    for leftover in re.findall(r"\[[^\]]*\]|\([^)]*\)", s):
        warnings.warn(
            f"Unrecognized modification {leftover!r} in sequence "
            f"{raw_seq!r}; skipping peptide.",
            stacklevel=2,
        )
        return None

    s = re.sub(r"\[.*?\]", "", s)
    s = re.sub(r"\(.*?\)", "", s)

    clean_chars = []
    oxM_positions = []

    for ch in s:
        if ch == "m":
            oxM_positions.append(len(clean_chars))
            clean_chars.append("M")
        elif ch.isalpha():
            clean_chars.append(ch.upper())

    clean_seq = "".join(clean_chars)
    oxM_count = len(oxM_positions)

    return clean_seq, oxM_count, oxM_positions

# Feature computation 
def compute_features(clean_seq, charge, oxM_count, oxM_positions, tmt=True):
    """Compute the 9 model features from sequence and charge"""

    L = len(clean_seq)
    z = float(charge)

    neutral_mass = pmass.calculate_mass(sequence=clean_seq)
    n_tmt = (1 + clean_seq.count("K")) if tmt else 0
    n_cys = clean_seq.count("C")
    precursor_mz = (neutral_mass
                    + n_tmt * TMTPRO_MASS
                    + oxM_count * OXM_MASS
                    + n_cys * CYS_MASS
                    + z * PROTON) / z

    frac_H = clean_seq.count("H") / L
    RK_frac = (clean_seq.count("R") + clean_seq.count("K")) / L

    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    pa = ProteinAnalysis(clean_seq)
    pI = pa.isoelectric_point()
    GRAVY = pa.gravy()

    if oxM_count > 0 and L > 1:
        oxM_position_mean = float(np.mean([p / (L - 1) for p in oxM_positions]))
    else:
        oxM_position_mean = 0.0

    pro_positions = [i for i, aa in enumerate(clean_seq) if aa == "P"]
    if pro_positions and L > 1:
        proline_position_mean = float(np.mean(pro_positions) / (L - 1))
    else:
        proline_position_mean = 0.0

    return {
        "charge": z,
        "precursor_mz": precursor_mz,
        "frac_H": frac_H,
        "oxM_count": float(oxM_count),
        "oxM_position_mean": oxM_position_mean,
        "RK_frac": RK_frac,
        "pI": pI,
        "GRAVY": GRAVY,
        "proline_position_mean": proline_position_mean,
    }

# Model loading & prediction 
def load_ensemble(model_dir, model_tag=MODEL_TAG):
    """Load the 5-fold XGBoost ensemble."""
    models = []
    for i in range(N_FOLDS):
        path = Path(model_dir) / f"{model_tag}_fold{i}.json"
        if not path.exists():
            raise FileNotFoundError(f"Model file not found: {path}")
        clf = xgb.XGBClassifier()
        clf.load_model(str(path))
        models.append(clf)
    return models

def predict_ensemble(models, features):
    """Average predicted probabilities across the 5 fold models."""
    X = np.array([[features[f] for f in FEATURE_ORDER]])
    proba_sum = np.zeros(len(ALL_SETTINGS))
    for clf in models:
        proba_sum += clf.predict_proba(X)[0]
    return proba_sum / len(models)


def main():
    """Predict the optimal fragmentation setting for peptides."""
    parser = argparse.ArgumentParser(
        description="Predict the optimal fragmentation setting for peptides."
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input CSV with columns: sequence, charge (optional: peptide_id)"
    )
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=Path(__file__).parent.parent / "model",
        help="Directory containing the 5-fold model files."
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output CSV file."
    )

    args = parser.parse_args()
    df = pd.read_csv(args.input)

    for col in ("sequence", "charge"):
        if col not in df.columns:
            raise ValueError(
                f"Input CSV is missing required column: '{col}'"
            )
    models = load_ensemble(args.model_dir)

    rows = []
    n_skipped = 0

    for _, row in df.iterrows():

        parsed = parse_sequence(str(row["sequence"]))

        # Skip peptides with unrecognized mod.
        if parsed is None:
            n_skipped += 1
            continue

        clean_seq, oxM_count, oxM_positions = parsed

        features = compute_features(
            clean_seq,
            row["charge"],
            oxM_count,
            oxM_positions
        )

        proba = predict_ensemble(models, features)
        order = np.argsort(proba)[::-1]

        out_row = {
            "peptide_id": (
                row["peptide_id"]
                if "peptide_id" in df.columns
                else ""
            ),
            "sequence": row["sequence"],
            "charge": row["charge"],
            "clean_sequence": clean_seq,
            "predicted_best_setting": ALL_SETTINGS[order[0]],
            "predicted_best_prob": float(proba[order[0]]),
        }

        out_row.update({
            f"P_{s}": float(p)
            for s, p in zip(ALL_SETTINGS, proba)
        })

        rows.append(out_row)

    out_df = pd.DataFrame(rows)

    if "peptide_id" not in df.columns:
        out_df = out_df.drop(columns=["peptide_id"])

    out_df.to_csv(args.output, index=False)
    print(f"Wrote {len(out_df)} predictions to {args.output}")

    if n_skipped:
        print(
            f"Skipped {n_skipped} row(s) with unrecognized "
            f"modification notation (see warnings above)."
        )


if __name__ == "__main__":
    main()
