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
