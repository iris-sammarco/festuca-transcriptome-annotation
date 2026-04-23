#!/usr/bin/env python3
"""
trinotate_09_merge_annotations.py
Merge Trinotate, EggNOG, and Rfam annotation sources into one extended TSV.
"""

import pandas as pd
import re
import os

############################################
# FILES (EDIT THESE)
############################################
BASE_DIR = "/path/to/project/assembly/final_annotation"

trinotate_file = os.path.join(BASE_DIR, "Festuca_rubra.annotation_report_strict.xls")
eggnog_file    = os.path.join(BASE_DIR, "eggnog.emapper.annotations")
rfam_file      = os.path.join(BASE_DIR, "rfam.tblout")
output_file    = os.path.join(BASE_DIR, "Festuca_rubra_annotation_report.extended.tsv")

############################################
# 1. LOAD TRINOTATE
############################################
trin = pd.read_csv(trinotate_file, sep="\t")

############################################
# 2. PROCESS EGGNOG
############################################
print("Processing eggNOG...")

egg_header = [
    "query",
    "seed_ortholog",
    "evalue",
    "score",
    "eggNOG_OGs",
    "max_annot_lvl",
    "COG_category",
    "Description",
    "Preferred_name",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_Reaction",
    "PFAMs"
]

egg = pd.read_csv(
    eggnog_file,
    sep="\t",
    names=egg_header,
    skiprows=1, # skip the header line that starts with #
    low_memory=False
)

egg = egg.rename(columns={"query": "query_id"})

# map protein → transcript
egg["transcript_id"] = egg["query_id"].str.replace(r"\.p\d+$", "", regex=True)

# keep useful columns
egg = egg[[
    "transcript_id",
    "evalue",
    "Preferred_name",
    "Description",
    "GOs",
    "KEGG_ko",
    "KEGG_Pathway",
    "COG_category"
]]

# convert evalue
egg["evalue"] = pd.to_numeric(egg["evalue"], errors="coerce")

############################################
# ⭐ FILTERING (EGGNOG)
############################################

# 1. remove weak hits
egg = egg[egg["evalue"] < 1e-10]

# 2. remove completely uninformative rows
egg = egg[
    egg["Description"].notna() &
    (egg["Description"] != "-")
]

############################################
# keep best hit per transcript
############################################
egg = egg.sort_values("evalue").drop_duplicates("transcript_id")

# rename columns
egg = egg.rename(columns={
    "Preferred_name": "eggnog_preferred_name",
    "Description": "eggnog_description",
    "GOs": "eggnog_GO",
    "KEGG_ko": "eggnog_KEGG_ko",
    "KEGG_Pathway": "eggnog_KEGG_pathway",
    "COG_category": "eggnog_COG",
    "evalue": "eggnog_evalue"
})

############################################
# 3. PROCESS RFAM
############################################
print("Processing Rfam...")

rfam_data = []

with open(rfam_file) as f:
    for line in f:
        if line.startswith("#"):
            continue

        parts = line.strip().split(maxsplit=26)
        # maxsplit=26 ensures last field = full description

        if len(parts) < 27:
            continue  # skip malformed lines safely

        rfam_data.append({
            "transcript_id": parts[3],
            "rfam_family": parts[1],
            "rfam_description": parts[26],
            "rfam_evalue": parts[17],
            "score": parts[16]
        })

rfam = pd.DataFrame(rfam_data)

# convert types
rfam["rfam_evalue"] = pd.to_numeric(rfam["rfam_evalue"], errors="coerce")
rfam["score"] = pd.to_numeric(rfam["score"], errors="coerce")

############################################
# ⭐ NO EXTRA FILTERING (already used --cut_ga)
############################################

############################################
# BEST HIT
############################################
rfam_best = rfam.sort_values("rfam_evalue").drop_duplicates("transcript_id")

rfam_best = rfam_best.rename(columns={
    "rfam_family": "rfam_best_family",
    "rfam_description": "rfam_best_description",
    "rfam_evalue": "rfam_best_evalue"
})

############################################
# ALL HITS (concatenated)
############################################
rfam_all = rfam.copy()

rfam_all["combined"] = (
    rfam_all["rfam_family"] + "|" +
    rfam_all["rfam_description"] + "|" +
    rfam_all["rfam_evalue"].astype(str)
)

rfam_all = rfam_all.groupby("transcript_id")["combined"] \
    .apply(lambda x: ";".join(x)) \
    .reset_index()

rfam_all = rfam_all.rename(columns={
    "combined": "rfam_all_hits"
})

############################################
# 4. MERGE
############################################
print("Merging...")

merged = trin.merge(egg, how="left", on="transcript_id")
merged = merged.merge(rfam_best, how="left", on="transcript_id")
merged = merged.merge(rfam_all, how="left", on="transcript_id")

############################################
# 5. SAVE
############################################
merged.to_csv(output_file, sep="\t", index=False)

print(f"Done! Output written to: {output_file}")
