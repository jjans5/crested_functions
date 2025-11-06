import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def rowwise_corr(A: pd.DataFrame, B: pd.DataFrame, method: str = "pearson") -> pd.DataFrame:
    """
    Row-wise correlation between A (n x c) and B (n x c) across columns.
    Returns an (n x n) DataFrame R where R[i, j] = corr(A[i, :], B[j, :]).
    
    method: "pearson" (default) or "spearman"
    Notes:
      - Assumes no NaNs; for NaNs use a pairwise handler (slower).
      - Uses unbiased ddof=1, so we divide by (m-1) in the final dot.
    """
    # 1) Align on common columns (order is preserved from A)
    cols = A.columns.intersection(B.columns)
    if len(cols) == 0:
        raise ValueError("A and B have no columns in common.")
    A_ = A[cols]
    B_ = B[cols]
    
    # 2) Optional: Spearman ranks per row (rank across columns)
    if method.lower() == "spearman":
        A_ = A_.rank(axis=1, method="average")
        B_ = B_.rank(axis=1, method="average")
    elif method.lower() != "pearson":
        raise ValueError("method must be 'pearson' or 'spearman'")
    
    # 3) Convert to numpy
    X = A_.to_numpy(dtype=float)
    Y = B_.to_numpy(dtype=float)
    m = X.shape[1]  # number of columns used
    
    # 4) Center rows
    X0 = X - X.mean(axis=1, keepdims=True)
    Y0 = Y - Y.mean(axis=1, keepdims=True)
    
    # 5) Std per row (ddof=1); protect against zeros
    Xstd = X0.std(axis=1, ddof=1, keepdims=True)
    Ystd = Y0.std(axis=1, ddof=1, keepdims=True)
    Xstd[Xstd == 0] = 1.0
    Ystd[Ystd == 0] = 1.0
    
    # 6) Z-score rows
    Xz = X0 / Xstd
    Yz = Y0 / Ystd
    
    # 7) Correlation matrix via matrix multiply
    # For ddof=1 standardization, sum(z^2) = m-1 per row, so divide dot by (m-1)
    R = (Xz @ Yz.T) / (m - 1)
    
    # 8) Wrap as DataFrame
    R_df = pd.DataFrame(R, index=A.index, columns=B.index)
    return R_df

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def rowwise_corr(A: pd.DataFrame, B: pd.DataFrame, method: str = "pearson") -> pd.DataFrame:
    """
    Row-wise correlation between A (n x c) and B (n x c) across columns.
    Returns an (n x n) DataFrame R where R[i, j] = corr(A[i, :], B[j, :]).
    
    method: "pearson" (default) or "spearman"
    Notes:
      - Assumes no NaNs; for NaNs use a pairwise handler (slower).
      - Uses unbiased ddof=1, so we divide by (m-1) in the final dot.
    """
    # 1) Align on common columns (order is preserved from A)
    cols = A.columns.intersection(B.columns)
    if len(cols) == 0:
        raise ValueError("A and B have no columns in common.")
    A_ = A[cols]
    B_ = B[cols]
    
    # 2) Optional: Spearman ranks per row (rank across columns)
    if method.lower() == "spearman":
        A_ = A_.rank(axis=1, method="average")
        B_ = B_.rank(axis=1, method="average")
    elif method.lower() != "pearson":
        raise ValueError("method must be 'pearson' or 'spearman'")
    
    # 3) Convert to numpy
    X = A_.to_numpy(dtype=float)
    Y = B_.to_numpy(dtype=float)
    m = X.shape[1]  # number of columns used
    
    # 4) Center rows
    X0 = X - X.mean(axis=1, keepdims=True)
    Y0 = Y - Y.mean(axis=1, keepdims=True)
    
    # 5) Std per row (ddof=1); protect against zeros
    Xstd = X0.std(axis=1, ddof=1, keepdims=True)
    Ystd = Y0.std(axis=1, ddof=1, keepdims=True)
    Xstd[Xstd == 0] = 1.0
    Ystd[Ystd == 0] = 1.0
    
    # 6) Z-score rows
    Xz = X0 / Xstd
    Yz = Y0 / Ystd
    
    # 7) Correlation matrix via matrix multiply
    # For ddof=1 standardization, sum(z^2) = m-1 per row, so divide dot by (m-1)
    R = (Xz @ Yz.T) / (m - 1)
    
    # 8) Wrap as DataFrame
    R_df = pd.DataFrame(R, index=A.index, columns=B.index)
    return R_df
