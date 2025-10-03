# sa_score_combined.py
import sys
import pandas as pd
import numpy as np
from rdkit import Chem
from sascorer import calculateScore  # SA score
from standalone_model_numpy import SCScorer  # SCS score

def get_input_smiles(input_arg):
    """
    Parse input: can be a CSV file, a single SMILES, or a list of SMILES.
    Returns a pandas DataFrame with a 'SMILES' column.
    """
    # If ends with .csv, read as CSV
    if input_arg.lower().endswith(".csv"):
        df = pd.read_csv(input_arg)
        if 'SMILES' not in df.columns:
            raise ValueError("CSV must have a 'SMILES' column.")
        return df
    else:
        # Otherwise treat as a single SMILES or comma-separated SMILES
        smiles_list = [s.strip() for s in input_arg.split(",")]
        return pd.DataFrame({'SMILES': smiles_list})


# --- Get input ---
input_arg = sys.argv[1]
output_file = sys.argv[2]
df = get_input_smiles(input_arg)



# --- SA scores ---
sa_scores = []
invalid_sa = []

for smi in df['SMILES']:
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        sa = calculateScore(mol)
    else:
        sa = np.nan
        invalid_sa.append(smi)
    sa_scores.append(sa)

df['SA_score'] = sa_scores

# --- SCS scores ---
model = SCScorer()
model.restore("model.ckpt-10654.as_numpy.json.gz")

scs_scores = []

for smi in df['SMILES']:
    try:
        _, scs = model.get_score_from_smi(smi)
    except:
        scs = np.nan
    scs_scores.append(scs)

df['SCS_score'] = scs_scores

# --- Molecular complexity scores ---
def is_valid_smiles(smiles):
    return Chem.MolFromSmiles(smiles) is not None

def mc1(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.nan
    frac = sum(1 for a in mol.GetAtoms() if a.GetDegree() == 2) / mol.GetNumAtoms()
    return 1 - frac

def mc2(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.nan
    atoms_in_COX = set()
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            b, e = bond.GetBeginAtom(), bond.GetEndAtom()
            if {b.GetAtomicNum(), e.GetAtomicNum()} == {6, 8}:  # C=O
                carbon = b if b.GetAtomicNum() == 6 else e
                oxygen = e if carbon == b else b
                for n in carbon.GetNeighbors():
                    if n.GetIdx() != oxygen.GetIdx() and n.GetAtomicNum() in [7, 8]:
                        atoms_in_COX.update([carbon.GetIdx(), oxygen.GetIdx()])
                        break
    count = sum(1 for a in mol.GetAtoms() if a.GetDegree() != 2 and a.GetIdx() not in atoms_in_COX)
    return count

mc1_scores = []
mc2_scores = []
invalid_mc = []

for smi in df['SMILES']:
    if not is_valid_smiles(smi):
        invalid_mc.append(smi)
        mc1_scores.append(np.nan)
        mc2_scores.append(np.nan)
        continue
    try:
        mc1_scores.append(mc1(smi))
        mc2_scores.append(mc2(smi))
    except:
        mc1_scores.append(np.nan)
        mc2_scores.append(np.nan)

df['MC1_score'] = mc1_scores
df['MC2_score'] = mc2_scores

# --- Save combined results ---
df.to_csv(output_file, index=False)
print(f"All scores saved to {output_file}")

