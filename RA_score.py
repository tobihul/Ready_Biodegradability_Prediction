# RA_score_pipeline.py
import sys
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import pickle
import tensorflow as tf
from tqdm import tqdm

if len(sys.argv) < 3:
    print("Usage: python RA_score_pipeline.py input.csv output.csv")
    sys.exit(1)


input_file = sys.argv[1]
output_file = sys.argv[2]

# --- Read input CSV ---
df = pd.read_csv(input_file)

# --- NN scorer ---
class RAScorerNN:
    def __init__(self, model_path="RAscore/RAscore/models/DNN_chembl_fcfp_counts/model.tf"):
        self.nn_model = tf.keras.models.load_model(model_path)

    def ecfp_counts(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES")
        fp = AllChem.GetMorganFingerprint(mol, 3, useCounts=True, useFeatures=True)
        arr = np.zeros((2048,), dtype=np.int32)
        for idx, v in fp.GetNonzeroElements().items():
            arr[idx % 2048] += int(v)
        return arr

    def predict(self, smiles):
        try:
            arr = self.ecfp_counts(smiles)
            return float(self.nn_model.predict(arr.reshape(1, -1), verbose=0)[0][0])
        except:
            return np.nan

# --- XGB scorer ---
class RAScorerXGB:
    def __init__(self, model_path="RAscore/RAscore/models/XGB_chembl_ecfp_counts/model.pkl"):
        with open(model_path, "rb") as f:
            self.xgb_model = pickle.load(f)

    def ecfp(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES")
        fp = AllChem.GetMorganFingerprint(mol, 3, useCounts=True, useFeatures=False)
        arr = np.zeros((2048,), dtype=np.int32)
        for idx, v in fp.GetNonzeroElements().items():
            arr[idx % 2048] += int(v)
        return arr

    def predict(self, smiles):
        try:
            arr = self.ecfp(smiles)
            proba = self.xgb_model.predict_proba(arr.reshape(1, -1))
            return float(proba[0][1])
        except:
            return np.nan

# --- Initialize scorers ---
nn_scorer = RAScorerNN()
xgb_scorer = RAScorerXGB()

# --- Predict RA scores ---
nn_scores = []
xgb_scores = []

for smi in tqdm(df['SMILES'], desc="Calculating RA scores"):
    nn_scores.append(nn_scorer.predict(smi))
    xgb_scores.append(xgb_scorer.predict(smi))

df['RA_NN_score'] = nn_scores
df['RA_XGB_score'] = xgb_scores

# --- Save output ---
df.to_csv(output_file, index=False)
print(f"RA Scores saved to {output_file}")


