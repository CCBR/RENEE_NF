import subprocess


def test_help():
    output = subprocess.run(
        "bin/renee_nf --help", capture_output=True, shell=True, text=True
    ).stdout
    assert "RENEE_NF" in output


def test_version():
    output = subprocess.run(
        "bin/renee_nf --version", capture_output=True, shell=True, text=True
    ).stdout
    assert "renee_nf, version " in output


def test_citation():
    output = subprocess.run(
        "bin/renee_nf --citation", capture_output=True, shell=True, text=True
    ).stdout
    assert "@misc{" in output


def test_subcommands_help():
    assert all(
        [
            f"renee_nf {cmd} [OPTIONS]"
            in subprocess.run(
                f"bin/renee_nf {cmd} --help",
                capture_output=True,
                shell=True,
                text=True,
            ).stdout
            for cmd in ["run", "init"]
        ]
    )
