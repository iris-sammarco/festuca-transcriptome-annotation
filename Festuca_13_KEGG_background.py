#!/usr/bin/env python3

"""
Output file:
transcript_id    KEGG_all
TRINITY_...      K00001;K01810
"""
import pandas as pd
import re

INPUT_FILE = "/path/to/project/assembly/final_annotation/Festuca_rubra_annotation_report.extended.tsv"
df = pd.read_csv(INPUT_FILE, sep="\t")

def extract_kegg(text):
    if pd.isna(text):
        return []
    return re.findall(r"K\d{5}", str(text))

def merge_kegg(row):
    kegg = set()
    for col in ["Kegg", "eggnog_KEGG_ko"]:
        kegg.update(extract_kegg(row.get(col)))
    return ";".join(sorted(kegg))

df["KEGG_all"] = df.apply(merge_kegg, axis=1)

kegg_df = df[["transcript_id", "KEGG_all"]]
kegg_df = kegg_df[kegg_df["KEGG_all"] != ""]

kegg_df.to_csv("Festuca_rubra_kegg_background.tsv", sep="\t", index=False)

#######################################################
# Cleanup: Remove duplicate rows and collapse isoforms
#######################################################

df2 = pd.read_csv("Festuca_rubra_kegg_background.tsv", sep="\t")

#################################
# 1. Remove exact duplicates
#################################
df2 = df2.drop_duplicates()

#################################
# 2. Collapse isoforms → gene level (OPTIONAL but recommended)
#################################
df2["gene_id"] = df2["transcript_id"].str.replace(r"_i\d+$", "", regex=True)

df2 = (
    df2.groupby("gene_id")["KEGG_all"]
    .apply(lambda x: ";".join(sorted(set(";".join(x).split(";")))))
    .reset_index()
)

#################################
# 3. Remove empty entries
#################################
df2 = df2[df2["KEGG_all"] != ""]

#################################
# 4. Save
#################################
df2.to_csv("Festuca_rubra_kegg_background.cleaned.tsv", sep="\t", index=False)
