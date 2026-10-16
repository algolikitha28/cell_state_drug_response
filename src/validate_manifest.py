import csv
import re
import tarfile
from collections import Counter
from pathlib import Path

METADATA_PATH = Path("data/raw/geo_sample_metadata.csv")
ARCHIVE_PATH = Path("data/raw/GSE140440_RAW.tar")
OUTPUT_PATH = Path("data/processed/validated_sample_manifest.csv")


def main():
    if not METADATA_PATH.exists():
        raise FileNotFoundError(f"Missing metadata: {METADATA_PATH}")

    if not ARCHIVE_PATH.exists():
        raise FileNotFoundError(f"Missing archive: {ARCHIVE_PATH}")

    # Read the GEO metadata.
    with METADATA_PATH.open("r", newline="", encoding="utf-8-sig") as handle:
        metadata = list(csv.DictReader(handle))

    metadata_accessions = [row["geo_accession"] for row in metadata]
    metadata_counts = Counter(metadata_accessions)

    # Inspect the archive without extracting its files.
    with tarfile.open(ARCHIVE_PATH, "r") as archive:
        members = [
            member for member in archive.getmembers()
            if member.isfile()
            and member.name.lower().endswith(".counts.txt.gz")
        ]

    archive_rows = []
    archive_counts = Counter()

    for member in members:
        filename = Path(member.name).name
        match = re.match(r"(GSM\d+)_", filename)

        if not match:
            archive_rows.append({
                "filename": filename,
                "accession": "",
                "archive_size_bytes": member.size,
            })
            continue

        accession = match.group(1)
        archive_counts[accession] += 1
        archive_rows.append({
            "filename": filename,
            "accession": accession,
            "archive_size_bytes": member.size,
        })

    metadata_set = set(metadata_accessions)
    archive_set = set(archive_counts)

    missing_files = sorted(metadata_set - archive_set)
    missing_metadata = sorted(archive_set - metadata_set)
    duplicate_metadata = sorted(
        accession for accession, count in metadata_counts.items()
        if count > 1
    )
    duplicate_archive = sorted(
        accession for accession, count in archive_counts.items()
        if count > 1
    )

    # Match each metadata row to its archive filename.
    archive_lookup = {}
    for row in archive_rows:
        accession = row["accession"]
        if accession and accession not in archive_lookup:
            archive_lookup[accession] = row

    manifest = []
    for row in metadata:
        accession = row["geo_accession"]
        archive_row = archive_lookup.get(accession)

        result = dict(row)
        result["archive_filename"] = (
            archive_row["filename"] if archive_row else ""
        )
        result["archive_size_bytes"] = (
            archive_row["archive_size_bytes"] if archive_row else ""
        )
        result["file_match_status"] = (
            "MATCHED" if archive_counts[accession] == 1
            else "MISSING_FILE" if archive_counts[accession] == 0
            else "DUPLICATE_ARCHIVE_ACCESSION"
        )
        manifest.append(result)

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(manifest[0].keys()) if manifest else []

    with OUTPUT_PATH.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(manifest)

    print("\n=== ARCHIVE-METADATA VALIDATION ===")
    print(f"Metadata records: {len(metadata)}")
    print(f"Count files in archive: {len(members)}")
    print(f"Unique metadata accessions: {len(metadata_set)}")
    print(f"Unique archive accessions: {len(archive_set)}")
    print(f"Metadata accessions without files: {len(missing_files)}")
    print(f"Archive accessions without metadata: {len(missing_metadata)}")
    print(f"Duplicate metadata accessions: {len(duplicate_metadata)}")
    print(f"Duplicate archive accessions: {len(duplicate_archive)}")

    if missing_files:
        print("Missing count files for:", ", ".join(missing_files[:20]))
    if missing_metadata:
        print("Files without metadata:", ", ".join(missing_metadata[:20]))
    if duplicate_metadata:
        print("Duplicate metadata:", ", ".join(duplicate_metadata[:20]))
    if duplicate_archive:
        print("Duplicate archive accessions:", ", ".join(duplicate_archive[:20]))

    print("\n=== CELL LINE x RESPONSE ===")
    combinations = Counter(
        (row["cell_line"], row["docetaxel_response"])
        for row in metadata
    )
    for cell_line in sorted({row["cell_line"] for row in metadata}):
        for response in ("Sensitive", "Resistant"):
            print(f"{cell_line} | {response}: {combinations[(cell_line, response)]}")

    matched_count = sum(
        row["file_match_status"] == "MATCHED" for row in manifest
    )
    passed = (
        len(metadata) == len(members)
        and len(missing_files) == 0
        and len(missing_metadata) == 0
        and len(duplicate_metadata) == 0
        and len(duplicate_archive) == 0
        and matched_count == len(metadata)
    )

    print(f"\nManifest saved to: {OUTPUT_PATH.resolve()}")
    print(f"Overall file-matching check: {'PASS' if passed else 'CHECK REQUIRED'}")


if __name__ == "__main__":
    main()
