#!/usr/bin/env python3
"""
trinotate_10_go_background.py
Write a GO background file (GO_all) combining GO terms from BLASTX, BLASTP, Pfam, and EggNOG.
Output columns: transcript_id  GO_all
"""

import pandas as pd
import re
import os

# FILES (EDIT THIS)
INPUT_FILE = "/path/to/project/assembly/final_annotation/Festuca_rubra_annotation_report.extended.tsv"

df = pd.read_csv(INPUT_FILE, sep="\t")

def extract_go(text):
    if pd.isna(text) or not text:
        return []
    return re.findall(r"GO:\d+", str(text))
