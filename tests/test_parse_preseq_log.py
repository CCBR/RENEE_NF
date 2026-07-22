import subprocess


def run_parser(tmp_path, contents):
    log = tmp_path / "preseq.log"
    log.write_text(contents)
    result = subprocess.run(
        ["bin/parse_preseq_log.py", str(log)],
        capture_output=True,
        check=True,
        text=True,
    )
    return result.stdout.strip()


def test_parse_preseq_log(tmp_path):
    output = run_parser(
        tmp_path,
        "TOTAL READS = 100\nDISTINCT READS = 80\n1\t60\n2\t10\n",
    )
    assert output == "0.8\t0.75\t6.0"


def test_parse_preseq_log_with_missing_metrics_returns_na(tmp_path):
    output = run_parser(tmp_path, "ERROR: reads unsorted\n")
    assert output == "NA\tNA\tNA"


def test_parse_preseq_log_with_zero_denominator_returns_na(tmp_path):
    output = run_parser(
        tmp_path,
        "TOTAL READS = 0\nDISTINCT READS = 0\n1\t0\n2\t0\n",
    )
    assert output == "NA\tNA\tNA"
