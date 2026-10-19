# FragWise

## Overview
FragWise predicts the MS2 fragmentation setting that is most likely to improve downstream quantification for a given precursor. It provides a potential framework for integration with real-time search platforms, where peptide sequences identified during data acquisition could be used to guide precursor-specific fragmentation in real time.

The model predicts probabilities for 10 conditions:

- CID 22, 25, 28, 31, 34
- HCD 22, 25, 28, 31, 34

🐱 The setting with the highest probability is reported as the predicted best setting for each precursor.

## Quick start
FragWise requires Python 3.12.

```bash
conda create -n fragwise python=3.12
conda activate fragwise
pip install -r requirements.txt
```

The trained 5-fold models are included in `model/`.

Run:

```bash
python scripts/predict.py \
    --input example/example_input.csv \
    --output example/example_output.csv
```

## Input
CSV with one row per peptide precursor:

| Column       | Description                                                 |
| ------------ | ----------------------------------------------------------- |
| `sequence`   | Peptide sequence. Oxidized Met: `M*`, `M(ox)`, or `M[+16]`. |
| `charge`     | Precursor charge, e.g. 2, 3, 4.                             |
| `peptide_id` | Optional.                                                     |

See `example/example_input.csv`.

## Output
* `predicted_best_setting` — predicted best fragmentation setting
* `predicted_best_prob` — predicted probability
* `P_<setting>` — probabilities for all 10 settings
* `clean_sequence` — sequence used for feature calculation

Possible settings:

```text
CID_22  CID_25  CID_28  CID_31  CID_34
HCD_22  HCD_25  HCD_28  HCD_31  HCD_34
```

See `example/example_output.csv`.

## Important

**TMTpro only.** FragWise was trained on TMTpro-labeled precursors. Non-TMT and label-free data are not supported.
Peptides with unrecognized modification notation are skipped with a warning.
For reproducibility, use the versions in `requirements.txt`. In particular:

```text
xgboost==3.2.0
```
Older XGBoost versions may produce different predictions.


## Model

- **Algorithm:** XGBoost, 5-fold ensemble.
- **Target:** Predict the best MS2 fragmentation setting for each peptide precursor.
- **Conditions:** CID/HCD at NCE 22, 25, 28, 31, or 34.
- **Features:** 9 features from sequence and charge. 
- **Prediction:** Probabilities from the 5 models are averaged, and the setting with the highest probability is reported.

## Files
```
fragwise/
├── README.md              
├── requirements.txt       
├── LICENSE
├── model/
│   └── fragwise_fold0-4.json              5-fold XGBoost ensemble
├── scripts/
│   └── predict.py         reads CSV -> computes features -> predicts -> writes CSV
├── example/
│   ├── example_input.csv   10 example peptides                      
│   └── example_output.csv  
└── docs/
    └── workflow.md       
```

## License

MIT — see [LICENSE](LICENSE).
