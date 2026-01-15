"""
Calculate Effective Sample Size (ESS) from BEAST log files.

This module provides utilities to read BEAST log files, calculate ESS using
autocorrelation analysis, and report on convergence thresholds.
"""

import numpy as np
import pandas as pd
from typing import Union, Tuple, List, Optional


def effective_sample_size(
    x: Union[np.ndarray, List[float]], max_lag: Optional[int] = None
) -> float:
    """
    Calculate effective sample size using autocorrelation.

    This matches the R effectiveSize function behavior:
    - Demeans the data
    - Calculates autocorrelation up to max_lag
    - Sums positive autocorrelations
    - ESS = N / (1 + 2 * sum)

    Parameters
    ----------
    x : array-like
        Time series data
    max_lag : int, optional
        Maximum lag for autocorrelation calculation. If None, uses
        min(N-1, 10000) where N is the length of x.

    Returns
    -------
    ess : float
        Effective sample size
    """
    x = np.asarray(x, dtype=float)
    N = len(x)

    if N < 2:
        return 0.0

    # Demean the data
    x_demeaned = x - np.mean(x)

    # Set max_lag if not provided
    if max_lag is None:
        max_lag = min(N - 1, 10000)
    else:
        max_lag = min(max_lag, N - 1)

    # Calculate autocorrelation function (matching R's acf with type="correlation")
    # Autocorrelation at lag k: r_k = Cov(X_t, X_{t+k}) / Var(X)
    variance = np.var(x_demeaned, ddof=0)
    if variance == 0:
        return float("inf")  # Constant series has infinite ESS

    # Calculate autocorrelations
    acf_lags = []
    for lag in range(1, max_lag + 1):
        # Optimized covariance calculation using slice
        # Cov(X_t, X_{t+k}) = mean((X_t - mu)(X_{t+k} - mu))
        covariance = np.sum(x_demeaned[:-lag] * x_demeaned[lag:]) / (N - lag)
        acf_val = covariance / variance
        if acf_val <= 0:
            break
        acf_lags.append(acf_val)

    # Calculate ESS
    pos_sum = sum(acf_lags)
    ess = N / (1 + 2 * pos_sum)

    return ess


def find_ess_threshold(
    data: Union[np.ndarray, List[float]],
    threshold: float = 200,
    max_lag: Optional[int] = None,
) -> Optional[int]:
    """
    Find the minimum number of samples needed to achieve ESS >= threshold.

    Parameters
    ----------
    data : array-like
        Logged sequential data
    threshold : float
        ESS threshold to achieve
    max_lag : int, optional
        Maximum lag for autocorrelation calculation

    Returns
    -------
    n_samples : int or None
        Minimum number of samples needed to achieve ESS >= threshold, or None if threshold not reached
    """
    data = np.asarray(data)
    N = len(data)

    # Check progressively larger sample sizes
    for n in range(100, N + 1, 100):
        ess = effective_sample_size(data[:n], max_lag)
        if ess >= threshold:
            # Refine search in smaller increments
            for n_refined in range(max(100, n - 100), n + 1, 10):
                if effective_sample_size(data[:n_refined], max_lag) >= threshold:
                    return n_refined
            return n

    # Check full dataset
    if effective_sample_size(data, max_lag) >= threshold:
        return N

    return None


def read_log_file(filename: str) -> Tuple[List[str], np.ndarray]:
    """
    Read BEAST log file and extract data.

    This function parses the standard BEAST .log format, ignoring comments.

    Parameters
    ----------
    filename : str
        Path to log file

    Returns
    -------
    header : list
        Column names
    data : numpy.ndarray
        Data matrix (samples x parameters)
    """
    header = None
    data_rows = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue

            if header is None:
                header = line.split("\t")
                continue

            if line:
                values = line.split("\t")
                try:
                    row = [float(v) for v in values]
                    data_rows.append(row)
                except ValueError:
                    continue

    if header is None:
        raise ValueError("No header found in log file")

    if len(data_rows) == 0:
        raise ValueError("No data rows found in log file")

    data = np.array(data_rows)

    if data.shape[1] != len(header):
        raise ValueError(
            f"Data columns ({data.shape[1]}) don't match header ({len(header)})"
        )

    return header, data


def analyze_ess(
    log_source: Union[str, pd.DataFrame],
    output_path: str,
    burnin: float = 0.1,
    check_threshold: bool = True,
) -> pd.DataFrame:
    """
    Calculate ESS for a BEAST log file or DataFrame and save results.

    Parameters
    ----------
    log_source : str or pd.DataFrame
        Path to the BEAST log file or a pandas DataFrame.
    output_path : str
        Path to save the CSV results table.
    burnin : float, optional
        Fraction of samples to discard as burn-in (0.0 to 1.0), default 0.1.
    check_threshold : bool, optional
        If True, prints a report on how many samples were needed to reach ESS > 200
        for 'posterior', 'prior', and 'likelihood' columns (if present aka intersection).

    Returns
    -------
    pd.DataFrame
        DataFrame containing the calculated ESS values for each parameter.
    """
    # 1. Load Data
    if isinstance(log_source, str):
        # BEAST logs are typically tab-separated with '#' comments
        try:
            df = pd.read_csv(log_source, sep="\t", comment="#")
        except Exception as e:
            raise ValueError(f"Failed to read log file {log_source}: {e}")
    elif isinstance(log_source, pd.DataFrame):
        df = log_source
    else:
        raise TypeError("log_source must be a file path (str) or pandas DataFrame")

    n_samples = len(df)
    if n_samples == 0:
        raise ValueError("Input data is empty")

    # 2. Apply Burn-in
    if not 0 <= burnin < 1:
        raise ValueError("Burn-in must be in range [0, 1)")

    burnin_idx = int(n_samples * burnin)
    df_burnin = df.iloc[burnin_idx:].reset_index(drop=True)

    print(
        f"Total samples: {n_samples}. After {burnin * 100:.1f}% burn-in: {len(df_burnin)} samples."
    )

    # 3. Calculate ESS for numeric columns
    # Exclude 'sample' or 'Sample' column if present, though effective_sample_size handles linear trend okay,
    # usually we skip it.
    skip_cols = {"sample", "Sample", "state", "State"}
    numeric_cols = df_burnin.select_dtypes(include=[np.number]).columns

    results = []
    for col in numeric_cols:
        if col in skip_cols:
            continue

        # Calculate ESS
        val = effective_sample_size(df_burnin[col].values)
        results.append({"Parameter": col, "ESS": val})

    results_df = pd.DataFrame(results).sort_values("ESS")

    # 4. Save Results
    results_df.to_csv(output_path, index=False)
    print(f"ESS results saved to {output_path}")

    # 5. Check Thresholds (Optional)
    if check_threshold:
        target_cols = set(["posterior", "prior", "likelihood"])
        # Case insensitive matching to be robust? The user used exact casing in their snippet.
        # But let's check exact match first as per user snippet.
        cols_to_check = target_cols.intersection(df.columns)

        if cols_to_check:
            print("\nFinding samples needed for ESS > 200:")
            threshold = 200.0
            for col in cols_to_check:
                data_series = df[col].values
                n_needed = find_ess_threshold(data_series, threshold=threshold)

                if n_needed is not None:
                    print(f"{col}: ESS > {threshold} achieved after {n_needed} samples")
                else:
                    final_ess = effective_sample_size(df_burnin[col].values)
                    print(
                        f"{col}: ESS > {threshold} not achieved (final ESS: {final_ess:.2f})"
                    )
        else:
            print(
                "\ncheck_threshold=True but 'posterior', 'prior', or 'likelihood' columns not found."
            )

    return results_df
