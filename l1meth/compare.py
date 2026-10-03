"""Head-to-head: bulk-average versus locus-resolved L1 methylation.

Tests the published negative result directly. A multi-cancer study of L1
promoter methylation at 11 loci reported that locus-specific markers showed no
significant clinical advantage over global markers for separating tumour from
normal tissue. That study measured few loci with a noisy assay; this asks under
what regime -- how many loci carry signal, at what read depth -- the conclusion
would or would not still hold.

Both feature sets are evaluated with the same estimator and the same folds, so
the comparison isolates the representation, not the classifier.
"""

from __future__ import annotations

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def _cv_auc(X: np.ndarray, y: np.ndarray, *, n_splits: int = 5, seed: int = 0) -> float:
    """Cross-validated ROC AUC for L2 logistic regression."""
    X = X.reshape(len(y), -1)
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    scores = np.empty(len(y), dtype=float)
    for tr, te in cv.split(X, y):
        model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
        model.fit(X[tr], y[tr])
        scores[te] = model.predict_proba(X[te])[:, 1]
    return float(roc_auc_score(y, scores))


def compare(X: np.ndarray, y: np.ndarray, *, seed: int = 0) -> dict[str, float]:
    """AUC for the bulk mean versus the full per-locus vector."""
    global_feature = X.mean(axis=1)
    return {
        "auc_global": _cv_auc(global_feature, y, seed=seed),
        "auc_locus": _cv_auc(X, y, seed=seed),
        "n_loci": X.shape[1],
        "n_samples": len(y),
    }
