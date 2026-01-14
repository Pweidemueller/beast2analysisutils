import numpy as np
from beast2analysisutils import effective_sample_size, read_log_file, analyze_ess


def test_effective_sample_size_constant():
    """Test ESS for constant data should be infinite (or handle 0 variance)."""
    # Our implementation returns inf for 0 variance
    x = np.ones(100)
    assert effective_sample_size(x) == float("inf")


def test_effective_sample_size_random():
    """Random uncorrelated data should have ESS close to N."""
    np.random.seed(42)
    N = 1000
    x = np.random.randn(N)
    ess = effective_sample_size(x)
    # For white noise, ESS should be close to N.
    # Let's check it's within 20% of N.
    assert N * 0.8 < ess < N * 1.2


def test_read_log_file_mock(tmp_path):
    """Test reading a mock BEAST log file."""
    log_file = tmp_path / "test.log"
    log_file.write_text(
        "#BEAST v2.5\nSample\tposterior\tlikelihood\n0\t-10\t-5\n1000\t-9\t-4"
    )

    header, data = read_log_file(str(log_file))
    assert header == ["Sample", "posterior", "likelihood"]
    assert data.shape == (2, 3)
    assert data[0, 1] == -10.0


def test_analyze_ess_wrapper(tmp_path):
    """Test the high-level analyze_ess wrapper."""
    log_file = tmp_path / "test.log"
    # Create enough data to not fail strict checks if any
    content = ["Sample\tposterior\tprior"]
    for i in range(101):
        content.append(f"{i * 1000}\t{np.random.randn()}\t{np.random.randn()}")
    log_file.write_text("\n".join(content))

    output_csv = tmp_path / "ess_results.csv"

    df_results = analyze_ess(
        str(log_file), str(output_csv), burnin=0.1, check_threshold=True
    )

    assert "Parameter" in df_results.columns
    assert "ESS" in df_results.columns
    assert "posterior" in df_results["Parameter"].values
    assert "prior" in df_results["Parameter"].values

    # Check file was created
    assert output_csv.exists()
