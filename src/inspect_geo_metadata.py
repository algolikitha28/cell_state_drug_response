import csv
import gzip
import re
from pathlib import Path
from collections import Counter

SOFT_PATH = Path("data/raw/GSE140440_family.soft.gz")
OUTPUT_PATH = Path("data/raw/geo_sample_metadata.csv")


def parse_geo_soft(path):
    samples = []
    current = None

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")

            # Start a new GEO sample record
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    samples.append(current)

                accession = line.split("=", 1)[1].strip()
                current = {
                    "geo_accession": accession,
                    "characteristics": [],
                    "data_processing": [],
                }
                continue

            if current is None or not line.startswith("!Sample_"):
                continue

            if " = " not in line:
                continue

            key, value = line.split(" = ", 1)
            value = value.strip()

            if key == "!Sample_title":
                current["title"] = value

            elif key == "!Sample_source_name_ch1":
                current["source_name"] = value

            elif key == "!Sample_characteristics_ch1":
                current["characteristics"].append(value)

            elif key == "!Sample_supplementary_file_1":
                current["supplementary_file"] = value

            elif key == "!Sample_treatment_protocol_ch1":
                current["treatment_protocol"] = value

            elif key == "!Sample_organism_ch1":
                current["organism"] = value

            elif key == "!Sample_data_processing":
                current["data_processing"].append(value)

    if current is not None:
        samples.append(current)

    return samples


def characteristic_value(characteristics, field_name):
    """Extract a value such as 'cell line: DU145'."""
    prefix = field_name.lower() + ":"

    for item in characteristics:
        key, separator, value = item.partition(":")
        if separator and key.strip().lower() == field_name.lower():
            return value.strip()

    return ""


def extract_sample(sample):
    characteristics = sample.get("characteristics", [])

    response = characteristic_value(characteristics, "docetaxel")
    cell_line = characteristic_value(characteristics, "cell line")
    disease = characteristic_value(characteristics, "disease state")

    # Use an explicit characteristic when available.
    if response.lower() not in {"sensitive", "resistant"}:
        response = "UNRESOLVED"

    supplementary_url = sample.get("supplementary_file", "")
    supplementary_filename = (
        supplementary_url.rsplit("/", 1)[-1]
        if supplementary_url else ""
    )

    # Extract the genome build only from an explicit processing statement.
    processing_text = " ".join(sample.get("data_processing", []))
    genome_match = re.search(
        r"Genome_build:\s*([^;,\s]+)",
        processing_text,
        flags=re.IGNORECASE,
    )

    return {
        "geo_accession": sample.get("geo_accession", ""),
        "title": sample.get("title", ""),
        "source_name": sample.get("source_name", ""),
        "cell_line": cell_line,
        "docetaxel_response": response,
        "disease_state": disease,
        "organism": sample.get("organism", ""),
        "supplementary_filename": supplementary_filename,
        "treatment_protocol": sample.get("treatment_protocol", ""),
        "genome_build": genome_match.group(1) if genome_match else "",
        "characteristics_all": " | ".join(characteristics),
    }


def main():
    if not SOFT_PATH.exists():
        raise FileNotFoundError(
            f"Metadata file not found: {SOFT_PATH}"
        )

    print("Parsing GEO family SOFT metadata...")

    raw_samples = parse_geo_soft(SOFT_PATH)
    records = [extract_sample(sample) for sample in raw_samples]

    fieldnames = [
        "geo_accession",
        "title",
        "source_name",
        "cell_line",
        "docetaxel_response",
        "disease_state",
        "organism",
        "supplementary_filename",
        "treatment_protocol",
        "genome_build",
        "characteristics_all",
    ]

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    with OUTPUT_PATH.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    print(f"Sample records parsed: {len(records)}")
    print(f"Metadata saved to: {OUTPUT_PATH.resolve()}")

    print("\nFirst 10 samples:")
    for row in records[:10]:
        print(
            f"{row['geo_accession']} | "
            f"{row['title']} | "
            f"{row['cell_line']} | "
            f"{row['docetaxel_response']} | "
            f"{row['supplementary_filename']}"
        )

    print("\nCell-line counts:")
    print(dict(Counter(r["cell_line"] or "MISSING" for r in records)))

    print("\nDocetaxel-response counts:")
    print(dict(Counter(r["docetaxel_response"] for r in records)))

    unresolved = [
        r["geo_accession"]
        for r in records
        if r["cell_line"] == ""
        or r["docetaxel_response"] == "UNRESOLVED"
    ]

    print(f"\nSamples with unresolved labels: {len(unresolved)}")
    if unresolved:
        print("Review these accessions:", ", ".join(unresolved[:20]))


if __name__ == "__main__":
    main()
