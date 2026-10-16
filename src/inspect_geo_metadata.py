from pathlib import Path
import csv
import gzip
import io

import requests


PROJECT_DIR = Path(__file__).resolve().parent.parent
RAW_DIR = PROJECT_DIR / "data" / "raw"

GEO_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/"
    "GSE140nnn/GSE140440/soft/"
    "GSE140440_family.soft.gz"
)

SOFT_PATH = RAW_DIR / "GSE140440_family.soft.gz"
OUTPUT_PATH = RAW_DIR / "geo_sample_metadata.csv"


def download_soft():
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    if SOFT_PATH.exists():
        print(f"Using existing file: {SOFT_PATH}")
        return

    print("Downloading GEO family SOFT metadata...")

    with requests.get(
        GEO_URL, stream=True, timeout=(30, 120)
    ) as response:
        response.raise_for_status()

        temp_path = SOFT_PATH.with_suffix(".gz.part")

        try:
            with temp_path.open("wb") as handle:
                for chunk in response.iter_content(1024 * 1024):
                    if chunk:
                        handle.write(chunk)

            # Check that the downloaded file is valid gzip.
            with gzip.open(temp_path, "rb") as handle:
                handle.read(1)

            temp_path.replace(SOFT_PATH)

        except Exception:
            temp_path.unlink(missing_ok=True)
            raise

    print("Download completed.")


def parse_soft():
    records = []
    current = None

    with gzip.open(
        SOFT_PATH, "rt", encoding="utf-8"
    ) as handle:
        for line in handle:
            line = line.rstrip("\r\n")

            if line.startswith("^SAMPLE = "):
                if current is not None:
                    records.append(current)

                current = {
                    "geo_accession": line.split("=", 1)[1].strip(),
                    "title": "",
                    "source_name_ch1": "",
                    "characteristics_ch1": [],
                }
                continue

            if current is None or not line.startswith("!Sample_"):
                continue

            # SOFT fields are tab-separated and can contain
            # quoted values. Parse them without losing tabs.
            parts = next(csv.reader([line], delimiter="\t"))
            key = parts[0]

            values = [
                value.strip()
                for value in parts[1:]
            ]

            if key == "!Sample_title":
                current["title"] = " | ".join(values)

            elif key == "!Sample_source_name_ch1":
                current["source_name_ch1"] = " | ".join(values)

            elif key.startswith("!Sample_characteristics_ch1"):
                current["characteristics_ch1"].extend(values)

        if current is not None:
            records.append(current)

    if not records:
        raise ValueError(
            "No sample records were found in the GEO SOFT file."
        )

    fieldnames = [
        "geo_accession",
        "title",
        "source_name_ch1",
        "characteristics_ch1",
    ]

    with OUTPUT_PATH.open(
        "w", newline="", encoding="utf-8-sig"
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames
        )
        writer.writeheader()

        for record in records:
            record["characteristics_ch1"] = " ; ".join(
                record["characteristics_ch1"]
            )
            writer.writerow(record)

    print(f"Sample records parsed: {len(records)}")
    print(f"Metadata saved to: {OUTPUT_PATH}")
    print("\nFirst 8 sample records:")

    for record in records[:8]:
        print(
            record["geo_accession"],
            "|", record["title"],
            "|", record["source_name_ch1"],
            "|", " ; ".join(record["characteristics_ch1"]),
        )


def main():
    download_soft()
    parse_soft()


if __name__ == "__main__":
    main()
