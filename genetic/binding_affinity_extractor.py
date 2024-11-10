import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.gaussian_process.kernels import Kernel
import joblib

class TanimotoKernel(Kernel):
    def __init__(self):
        pass

    def __call__(self, X, Y=None):
        if Y is None:
            Y = X
        K = np.zeros((X.shape[0], Y.shape[0]))
        for i in range(X.shape[0]):
            for j in range(Y.shape[0]):
                K[i, j] = self._tanimoto_similarity(X[i], Y[j])
        return K

    def _tanimoto_similarity(self, fp1, fp2):
        intersection = np.dot(fp1, fp2)
        union = np.sum(fp1) + np.sum(fp2) - intersection
        return intersection / union if union != 0 else 0.0

    def diag(self, X):
        return np.ones(X.shape[0])

    def is_stationary(self):
        return False

morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def smiles_to_fingerprint(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        mol = None

    if mol is not None:
        fp = morgan_gen.GetFingerprint(mol)
        return np.array(fp)
    return None

def load_novel_smiles(filepath):
    data = pd.read_csv(filepath, delimiter=';') 
    smiles_list = data['Smiles'].tolist()

    fingerprints = []
    valid_smiles = []
    invalid_smiles = []

    for smiles in smiles_list:
        fp = smiles_to_fingerprint(smiles)
        if fp is not None:
            valid_smiles.append(smiles)
            fingerprints.append(fp)
        else:
            invalid_smiles.append(smiles) 

    if invalid_smiles:
        print(f"Warning: Invalid SMILES found and ignored: {invalid_smiles}")

    fingerprints = np.array(fingerprints)
    return valid_smiles, fingerprints

def get_binding_affinities(smiles_filepath, model_filepath):
    gp = joblib.load(model_filepath)
    novel_smiles, X_novel = load_novel_smiles(smiles_filepath)

    if len(X_novel) == 0:
        raise ValueError("No valid fingerprints found.")

    predicted_affinities = gp.predict(X_novel)


    return predicted_affinities
