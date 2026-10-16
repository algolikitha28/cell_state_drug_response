from pathlib import Path
import gzip
import io
import tarfile

PROJECT_DIR = Path(__file__).resolve().parent.parent
ARCHIVE = PROJECT_DIR / "data" / "raw" / "GSE140440_RAW.tar"


def main():
    if not ARCHIVE.exists():
        raise FileNotFoundError(f"Archive not found: {ARCHIVE}")

    with tarfile.open(ARCHIVE, "r:") as tar:
        members = [
            m for m in tar.getmembers()
            if m.isfile() and m.name.endswith(".counts.txt.gz")
        ]

        print(f"Compressed count files: {len(members)}")
        print("\nFirst 10 filenames:")

        for member in members[:10]:
            print(f"  {member.name}")

        if not members:
            raise ValueError("No compressed count files found.")

        print("\nInspecting the first count file:")
        print(members[0].name)

        extracted = tar.extractfile(members[0])
        if extracted is None:
            raise OSError("Could not read the selected archive member.")

        with extracted:
            with gzip.GzipFile(fileobj=extracted) as compressed:
                with io.TextIOWrapper(
                    compressed, encoding="utf-8"
                ) as text_file:
                    for line_number, line in enumerate(text_file, 1):
                        print(f"{line_number}: {line.rstrip()}")
                        if line_number >= 10:
                            break


if __name__ == "__main__":
    main()
