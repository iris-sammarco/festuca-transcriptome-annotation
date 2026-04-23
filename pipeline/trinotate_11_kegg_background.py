#!/usr/bin/env python3
"""
trinotate_11_kegg_background.py
Write a KEGG background file combining KEGG terms from Trinotate and EggNOG.
Output columns: transcript_id  KEGG_all
"""

import pandas as pd
import re
import os

# FILES (EDIT THIS)
INPUT_FILE = "/path/to/project/assembly/final_annotation/Festuca_rubra_annotation_report.extended.tsv"

df = pd.read_csv(INPUT_FILE, sep="\t")

def extract_kegg(text):
    if pd.isna(text):
        return []
    return re.findall(r"K\d{5}", str(text))
