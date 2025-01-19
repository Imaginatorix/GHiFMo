import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.gaussian_process.kernels import Kernel
import joblib
from sklearn.preprocessing import MinMaxScaler
from Chromosome import *  

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
        if mol is not None:
            fp = morgan_gen.GetFingerprint(mol)
            return np.array(fp)
    except Exception:
        return None
    return None


def process_smiles(smiles_list):
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
        print(f"Ignored invalid SMILES: {invalid_smiles}")

    fingerprints = np.array(fingerprints)
    return valid_smiles, fingerprints


def get_binding_affinities(smiles_list, model_filepath, scaler_filepath):
    gp = joblib.load(model_filepath)  
    scaler = joblib.load(scaler_filepath)  

   
    valid_smiles, X_novel = process_smiles(smiles_list)

    if len(X_novel) == 0:
        raise ValueError("No valid fingerprints found.")

    X_novel_scaled = scaler.transform(X_novel)

    
    predicted_affinities_scaled = gp.predict(X_novel_scaled)

    predicted_affinities = scaler.inverse_transform(predicted_affinities_scaled.reshape(-1, 1)).flatten()

    return valid_smiles, predicted_affinities


def get_binding_affinity(molecules):
    smiles = [atom_to_smiles(molecule.head_atom) for molecule in molecules]
    model_filepath = "../genetic/trainer.pkl"  
    scaler_filepath = "../genetic/scaler.pkl"  

    valid_smiles, predicted_affinity = get_binding_affinities(smiles, model_filepath, scaler_filepath)
    return valid_smiles, predicted_affinity

