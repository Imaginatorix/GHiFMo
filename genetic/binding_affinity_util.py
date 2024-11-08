from sklearn.gaussian_process.kernels import Kernel
import numpy as np

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