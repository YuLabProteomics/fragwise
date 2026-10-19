# FragWise Workflow

This document provides the technical details needed to reproduce the FragWise prediction workflow.

## 1. System and software

FragWise was developed and tested with Python 3.12.

The package versions used for the manuscript are listed in `requirements.txt`:

```text
numpy==2.4.6
pandas==2.3.3
xgboost==3.2.0
biopython==1.86
pyteomics==4.7.5
scikit-learn==1.9.1
```

### System requirements

* **Linux:** No additional system packages are required.
* **macOS:** Install OpenMP first:

  ```bash
  brew install libomp
  ```
* **Windows:** Microsoft Visual C++ Redistributable is required.
* **GPU:** Not required.

`xgboost==3.2.0` is required for reproducibility. Older versions may load the released model but produce different predictions.

## 2. Sequence preprocessing

The script supports three oxidized-Met notations:

```text
M*
M(ox)
M[+16]
```

Before prediction, these annotations are removed from the sequence. The number and positions of oxidized methionines are kept as separate features. Unrecognized modifications are skipped with a warning. No normalization or scaling is applied.

## 3. Feature calculation

The model uses nine features:

| Feature                 | Description                                      |
| ----------------------- | ------------------------------------------------ |
| `charge`                | Precursor charge state                           |
| `precursor_mz`          | Theoretical precursor m/z                        |
| `frac_H`                | Fraction of histidine residues                   |
| `oxM_count`             | Number of oxidized methionines                   |
| `oxM_position_mean`     | Mean normalized position of oxidized methionines |
| `RK_frac`               | Combined fraction of arginine and lysine         |
| `pI`                    | Peptide isoelectric point                        |
| `GRAVY`                 | Grand Average of Hydropathicity                  |
| `proline_position_mean` | Mean normalized position of proline residues     |

### Precursor m/z

`precursor_mz` is calculated from the theoretical peptide mass using the modifications:

* TMTpro labeling at the peptide N-terminus and lysine residues
* carbamidomethylation of cysteine
* oxidation of methionine when annotated

TMTpro labeling is assumed for all predictions.

All nine features are calculated directly from the sequence and charge in `scripts/predict.py`.

## 4. Model details

FragWise uses five XGBoost models. All five models predict the same 10 fragmentation settings:

```text
0  CID_22
1  CID_25
2  CID_28
3  CID_31
4  CID_34
5  HCD_22
6  HCD_25
7  HCD_28
8  HCD_31
9  HCD_34
```

The XGBoost parameters are:

```text
n_estimators=200
max_depth=6
learning_rate=0.1
min_child_weight=10
subsample=0.8
colsample_bytree=0.8
reg_alpha=0.0
reg_lambda=1.0
```

The five models output probabilities for the 10 conditions. These probabilities are averaged, and the setting with the highest average is reported as `predicted_best_setting`.

## 5. Training data and scope

The released model was trained on TMTpro-labeled tryptic peptide precursors.

The training data included:

* charge states 2–5
* peptide lengths of 7–46 residues

The model only uses sequence- and charge-derived features. 

## 6. Reproducibility

The trained models are in `model/`, and the prediction script is in `scripts/predict.py`. Use the example files in `example/` to test the workflow. Use the package versions in `requirements.txt` for the same results.

## 7. Applying the model

To predict new peptides, provide a CSV with `sequence` and `charge` columns and run the prediction script as described in the README. The script calculates all nine features automatically. No retraining is needed.
