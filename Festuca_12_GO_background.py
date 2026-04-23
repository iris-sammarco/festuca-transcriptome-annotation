#!/usr/bin/env python3

"""
Writes a GO background file (GO_all) combining GO terms from BLASTX, BLASTP, Pfam, and eggNOG.
The file will look like:
transcript_id    GO_all
TRINITY_...      GO:0008150;GO:0003674;...
"""

import pandas as pd
import re

INPUT_FILE = "/path/to/project/assembly/final_annotation/Festuca_rubra_annotation_report.extended.tsv"
df = pd.read_csv(INPUT_FILE, sep="\t")

def extract_go(text):
    if pd.isna(text) or not text:
        return []
    return re.findall(r"GO:\d+", str(text))

def merge_go(row):
    gos = set()
    for col in [
        "gene_ontology_BLASTX",
        "gene_ontology_BLASTP",
        "gene_ontology_Pfam",
        "eggnog_GO"
    ]:
        gos.update(extract_go(row.get(col)))
    return ";".join(sorted(gos))

df["GO_all"] = df.apply(merge_go, axis=1)

# keep only rows with GO
go_df = df[["transcript_id", "GO_all"]]
go_df = go_df[go_df["GO_all"] != ""]

go_df.to_csv("Festuca_rubra_go_background.tsv", sep="\t", index=False)

############################################################################################
# Cleanup: Remove duplicate rows, collapse isoforms and remove generic GO terms (BP, MF, CC)
############################################################################################

df2 = pd.read_csv("Festuca_rubra_go_background.tsv", sep="\t")

#################################
# 1. Collapse isoforms → gene level
#################################
df2["gene_id"] = df2["transcript_id"].str.replace(r"_i\d+$", "", regex=True)

#################################
# 2. Merge GO terms per gene
#################################
def merge_go(series):
    gos = set()
    for entry in series.dropna():
        gos.update(entry.split(";"))
    return ";".join(sorted(gos))

df2 = (
    df2.groupby("gene_id")["GO_all"]
    .apply(merge_go)
    .reset_index()
)

#################################
# 3. Remove empty entries
#################################
df2 = df2[df2["GO_all"] != ""]

#################################
# 4. OPTIONAL: remove generic GO terms
#################################
GENERIC_GO = {
    "GO:0008150",  # biological_process
    "GO:0003674",  # molecular_function
    "GO:0005575"   # cellular_component
}

def filter_go(go_string):
    gos = [g for g in go_string.split(";") if g not in GENERIC_GO]
    return ";".join(gos)

df2["GO_all"] = df2["GO_all"].apply(filter_go)

#################################
# 5. Remove rows that became empty after filtering
#################################
df2 = df2[df2["GO_all"] != ""]

#################################
# 6. Save
#################################
df2.to_csv("Festuca_rubra_go_background.cleaned.tsv", sep="\t", index=False)
