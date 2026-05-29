import subprocess


def test_help():
    output = subprocess.run(
        "bin/renee_nf_alec --help", capture_output=True, shell=True, text=True
    ).stdout
    assert "RENEE_NF_ALEC" in output


def test_version():
    output = subprocess.run(
        "bin/renee_nf_alec --version", capture_output=True, shell=True, text=True
    ).stdout
    assert "renee_nf_alec, version " in output


def test_citation():
    output = subprocess.run(
        "bin/renee_nf_alec --citation", capture_output=True, shell=True, text=True
    ).stdout
    assert "@misc{" in output


def test_subcommands_help():
    assert all(
        [
            f"renee_nf_alec {cmd} [OPTIONS]"
            in subprocess.run(
                f"bin/renee_nf_alec {cmd} --help",
                capture_output=True,
                shell=True,
                text=True,
            ).stdout
            for cmd in ["run", "init"]
        ]
    )
