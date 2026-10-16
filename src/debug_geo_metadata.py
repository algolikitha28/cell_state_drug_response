import gzip

path = r"data/raw/GSE140440_family.soft.gz"
sample_to_find = "GSM4161452"
inside_sample = False

with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
    for line in handle:
        line = line.rstrip("\n")

        if line.startswith("^SAMPLE = "):
            accession = line.split("=", 1)[1].strip()
            inside_sample = accession == sample_to_find

            if inside_sample:
                print(f"\nSample block: {accession}")
                print("=" * 60)

            elif accession != sample_to_find:
                continue

        if inside_sample:
            if line.startswith("!Sample_") or line.startswith("^SAMPLE"):
                print(line)

            if line.startswith("!sample_table_begin"):
                break
