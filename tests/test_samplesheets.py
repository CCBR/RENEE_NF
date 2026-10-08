import csv
from pathlib import Path


def test_bundled_samplesheet_fastq_paths_resolve_from_project_directory():
    project_dir = Path(__file__).resolve().parents[1]
    samplesheets = (project_dir / "assets").glob("samplesheet*.csv")

    for samplesheet in samplesheets:
        with samplesheet.open(newline="") as csv_file:
            for row in csv.DictReader(csv_file):
                for fastq in (row["fastq_1"], row["fastq_2"]):
                    if fastq:
                        if Path(fastq).is_absolute():
                            continue
                        assert fastq.startswith("${projectDir}/")
                        fastq_path = Path(
                            fastq.replace("${projectDir}", str(project_dir))
                        )
                        assert fastq_path.resolve().is_file(), (
                            f"{fastq} from {samplesheet} does not resolve "
                            f"from {project_dir}"
                        )
