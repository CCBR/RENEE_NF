import csv
from pathlib import Path


def test_bundled_samplesheet_fastq_paths_resolve_from_init_directory():
    init_directory = Path("tmp")
    samplesheets = Path("assets").glob("samplesheet*.csv")

    for samplesheet in samplesheets:
        with samplesheet.open(newline="") as csv_file:
            for row in csv.DictReader(csv_file):
                for fastq in (row["fastq_1"], row["fastq_2"]):
                    if fastq:
                        fastq_path = Path(fastq)
                        if fastq_path.is_absolute():
                            continue
                        fastq_path = init_directory / fastq_path
                        assert fastq_path.resolve().is_file(), (
                            f"{fastq} from {samplesheet} does not resolve "
                            f"from {init_directory}"
                        )
