#to quantify differences between one reference sequence and 1/2 sequences
import argparse

from Bio import AlignIO
import pandas as pd

parser = argparse.ArgumentParser(
    description="Compare each aligned sequence against a reference."
)
parser.add_argument("alignment", help="Path to the aligned FASTA file")
parser.add_argument("reference_id", help="FASTA record ID of the reference")
parser.add_argument(
    "--out-prefix",
    default="",
    help="Prefix for both CSV files to distinguish between species/gene",
)
args = parser.parse_args()

alignment = AlignIO.read(args.alignment, "fasta")

seq_names = [r.id for r in alignment]
seqs = [list(str(r.seq)) for r in alignment]

#nseq = len(seqs)
#ncol = len(seqs[0])

# Running coordinates for each sequence
coords = []

for seq in seqs:

    running = 0
    coord_row = []

    for base in seq:

        if base != "-":
            running += 1
            coord_row.append(running)
        else:
            coord_row.append(None)

    coords.append(coord_row)

# Set this to the FASTA record ID of your reference.
reference_id = args.reference_id
if seq_names.count(reference_id) != 1:
    parser.error(f"Reference ID must appear exactly once in the FASTA: {reference_id!r}")
ref_idx = seq_names.index(reference_id)
reference = seqs[ref_idx]

rows = []
summary_rows = []

for query_idx, query in enumerate(seqs):
    if query_idx == ref_idx:
        continue

    query_id = seq_names[query_idx]
    shared = [
        col for col, (ref_base, query_base) in enumerate(zip(reference, query))
        if ref_base != "-" and query_base != "-"
    ]

    if not shared:
        print(f"{query_id}: no overlap with {reference_id}")
        continue

    mismatches = gaps = 0

    for col in range(shared[0], shared[-1] + 1):
        ref_base = reference[col].upper()
        query_base = query[col].upper()
        if ref_base == query_base:
            continue

        if "-" in (ref_base, query_base):
            difference = "gap"
            gaps += 1
        else:
            difference = "mismatch"
            mismatches += 1

        rows.append({
            "Reference": reference_id,
            "Query": query_id,
            "AlignPos": col + 1,
            "RefBase": ref_base,
            "RefPos": coords[ref_idx][col],
            "QueryBase": query_base,
            "QueryPos": coords[query_idx][col],
            "Difference": difference,
            "Type": "INDEL" if difference == "gap" else "substitution",
        })
    #  for a little summary printed out  
    compared_columns = sum(
        reference[col] != "-" or query[col] != "-"
        for col in range(shared[0], shared[-1] + 1)
    )
    total = mismatches + gaps

    summary_rows.append({
        "Reference": reference_id,
        "Query": query_id,
        "StartAlignPos": shared[0] + 1,
        "EndAlignPos": shared[-1] + 1,
        "ComparedColumns": compared_columns,
        "Substitutions": mismatches,
        "GapColumns": gaps,
        "TotalDifferences": total,
        "PercentDifferent": round(100 * total / compared_columns, 2),
    }) 
    #
    print(
        f"{query_id}: {mismatches} mismatches, {gaps} gap columns, "
        f"{mismatches + gaps} total differences"
    )

df = pd.DataFrame(rows, columns=[
    "Reference", "Query", "AlignPos", "RefBase", "RefPos",
    "QueryBase", "QueryPos", "Difference", "Type",
])
df = df.astype({"RefPos": "Int64", "QueryPos": "Int64"})

prefix = f"{args.out_prefix}_" if args.out_prefix else ""

df.to_csv(f"{prefix}variable_sites.csv", index=False)
pd.DataFrame(summary_rows).to_csv(f"{prefix}summary.csv", index=False)

print(df.head())

#
"""
if "-" in bases:
    variant_type = "INDEL"
else:
    variant_type = "SNP"

row["Type"] = variant_type
#
rows = []

for col in range(ncol):

    bases = [seq[col] for seq in seqs]

    if len(set(bases)) == 1:
        continue

    row = {
        "AlignPos": col + 1
    }

    for i, name in enumerate(seq_names):

        base = seqs[i][col]
        pos = coords[i][col]

        if base != "-":
            row[f"{name}_base"] = base
            row[f"{name}_pos"] = pos

        else:
            # nearest flanking position
            left = None

            for j in range(col - 1, -1, -1):
                if coords[i][j] is not None:
                    left = coords[i][j]
                    break

            if left is not None:
                row[f"{name}_base"] = "-"
                row[f"{name}_pos"] = f"{left}.5"
            else:
                row[f"{name}_base"] = "-"
                row[f"{name}_pos"] = "start_gap"

    rows.append(row)

df = pd.DataFrame(rows)

df.to_csv("variable_sites.csv", index=False)

print(df.head())
"""
