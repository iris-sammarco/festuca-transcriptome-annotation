#!/bin/bash

# parallel_trinotate_steps.bash - Split FASTA -> PBS array jobs -> merge outputs for Pfam/EggNOG
## NO QSUB for Master Script!
# Usage: bash trinotate_03_parallel_pfam_eggnog.bash pfam 141    # split+submit
#        bash trinotate_03_parallel_pfam_eggnog.bash pfam 141 m # merges the chunks, creates pfam.domtblout
#        bash trinotate_03_parallel_pfam_eggnog.bash eggnog 72 # split+submit
#        bash trinotate_03_parallel_pfam_eggnog.bash eggnog 72 m  # merges the chunks, creates eggnog.emapper.annotations

# =============================================================================
# trinotate_03_parallel_pfam_eggnog.bash - Parallel Trinotate Pfam/EggNOG Analysis
# =============================================================================
# AUTHOR: Iris Sammarco
#
# DATE: 25/03/2026
#
# PURPOSE:
#   Splits large TransDecoder peptide FASTA files into chunks and submits PBS array jobs for parallel Pfam domain search (hmmscan) or EggNOG annotation.
#   Includes merge step to combine chunked outputs into final annotation files.
#
# USAGE:
#   1. SPLIT + SUBMIT (creates chunks, submits PBS array):
#      bash trinotate_03_parallel_pfam_eggnog.bash pfam 141      # 141 Pfam chunks
#      bash trinotate_03_parallel_pfam_eggnog.bash eggnog 141    # 141 EggNOG chunks
#
#   2. MERGE ONLY (after jobs complete):
#      bash trinotate_03_parallel_pfam_eggnog.bash pfam 141 m    # Merge Pfam -> pfam.domtblout
#      bash trinotate_03_parallel_pfam_eggnog.bash eggnog 141 m  # Merge EggNOG -> eggnog.emapper.annotations
#
# ARRAY SCRIPTS CALLED:
#   trinotate_03_pfam_array.pbs
#   trinotate_03_eggnog_array.pbs

set -euo pipefail
trap 'echo "[ERROR] Line $LINENO: ${BASH_COMMAND} failed on $(date)" >&2; exit 1' ERR

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: $0 {pfam|eggnog} N_CHUNKS [m]" >&2
    exit 1
fi

export TMPDIR=$SCRATCHDIR
module add mambaforge
export CONDA_ENVS_PATH=/path/to/conda/envs
mamba activate trinotate_env

STEP=$1
N_CHUNKS=$2
MERGE_ONLY=${3:-}
OUTDIR="/path/to/project/assembly/trinotate_output"
PEP_LONG="${OUTDIR}/Trinity.fasta.transdecoder_dir/longest_orfs.pep"
PFAM_HMM="/path/to/project/assembly/trinotate_data/Pfam-A.hmm"

cd "${OUTDIR}"

case "${STEP}" in
    pfam)
        OUTPUT="pfam.domtblout"
        SPLIT_DIR="${OUTDIR}/Trinity.fasta.transdecoder_dir/longest_orfs.pep.split"
        PBS_SCRIPT="trinotate_03_pfam_array.pbs"
        INPUT_FILE="${PEP_LONG}"
        CHUNK_EXT="pep"
        ;;
    eggnog)
        OUTPUT="eggnog.emapper.annotations"
        SPLIT_DIR="${OUTDIR}/Trinity.fasta.transdecoder.pep.split"
        PBS_SCRIPT="trinotate_03_eggnog_array.pbs"
        INPUT_FILE="${OUTDIR}/Trinity.fasta.transdecoder.pep"
        CHUNK_EXT="pep"
        ;;
    *)
        echo "[ERROR] Unknown step: ${STEP}. Use pfam or eggnog" >&2
        exit 1
        ;;
esac
