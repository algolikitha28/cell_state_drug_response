from pathlib import Path
import csv
import hashlib
import tarfile

import requests


# Project paths
PROJECT_DIR = Path(__file__).resolve().parent.parent
RAW_DIR = PROJECT_DIR / "data" / "raw"

GEO_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/"
    "GSE140nnn/GSE140440/suppl/GSE140440_RAW.tar"
)

ARCHIVE_PATH = RAW_DIR / "GSE140440_RAW.tar"
CHECKSUM_PATH = RAW_DIR / "GSE140440_RAW.sha256"
MANIFEST_PATH = RAW_DIR / "archive_manifest.csv"


def calculate_sha256(file_path):
    """Calculate a file's SHA-256 checksum in small blocks."""
    digest = hashlib.sha256()

    with file_path.open("rb") as file:
        for block in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(block)

    return digest.hexdigest()


def download_archive():
    """Download the GEO archive without keeping incomplete downloads."""
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    if ARCHIVE_PATH.exists():
        print(f"Archive already exists: {ARCHIVE_PATH}")
        print("Calculating its checksum instead of downloading again.")
        return

    partial_path = ARCHIVE_PATH.with_suffix(".tar.part")

    try:
        print(f"Downloading from:\n{GEO_URL}")

        with requests.get(
            GEO_URL,
            stream=True,
            timeout=(30, 120),
        ) as response:
            response.raise_for_status()

            with partial_path.open("wb") as file:
                for chunk in response.iter_content(
                    chunk_size=1024 * 1024
                ):
                    if chunk:
                        file.write(chunk)

        # Rename only after the download finishes successfully.
        partial_path.replace(ARCHIVE_PATH)
        print("Download completed.")

    except Exception:
        partial_path.unlink(missing_ok=True)
        raise


def inspect_archive():
    """Validate the TAR and record its member names and sizes."""
    if not tarfile.is_tarfile(ARCHIVE_PATH):
        raise ValueError(
            "The downloaded file is not a valid TAR archive."
        )

    rows = []

    with tarfile.open(ARCHIVE_PATH, mode="r:") as archive:
        for member in archive.getmembers():
            # Record paths before extraction; do not extract yet.
            rows.append({
                "member_name": member.name,
                "size_bytes": member.size,
                "is_file": member.isfile(),
                "is_directory": member.isdir(),
            })

    with MANIFEST_PATH.open(
        "w", newline="", encoding="utf-8"
    ) as file:
        writer = csv.DictWriter(
            file,
            fieldnames=[
                "member_name",
                "size_bytes",
                "is_file",
                "is_directory",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Archive members: {len(rows)}")
    print(f"Manifest saved to: {MANIFEST_PATH}")

    # Print a small preview so we can inspect file naming.
    for row in rows[:15]:
        print(
            f"{row['member_name']} "
            f"({row['size_bytes']} bytes)"
        )


def save_checksum():
    """Save the archive's SHA-256 checksum."""
    checksum = calculate_sha256(ARCHIVE_PATH)

    CHECKSUM_PATH.write_text(
        f"{checksum}  {ARCHIVE_PATH.name}\n",
        encoding="utf-8",
    )

    print(f"SHA-256: {checksum}")
    print(f"Checksum saved to: {CHECKSUM_PATH}")


def main():
    download_archive()
    save_checksum()
    inspect_archive()


if __name__ == "__main__":
    main()
