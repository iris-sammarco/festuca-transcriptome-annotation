#!/bin/bash

# Run with:
# bash trinotate_05_parallel_signalp6.bash 72   # Submit chunks
# bash trinotate_05_parallel_signalp6.bash 72 m # Merge and load after all jobs finish

set -euo pipefail
trap 'echo "[ERROR] Line $LINENO: ${BASH_COMMAND} failed on $(date)" >&2; exit 1' ERR

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 N_CHUNKS [m]" >&2
    exit 1
fi

export TMPDIR=$SCRATCHDIR
module add mambaforge
export CONDA_ENVS_PATH=/path/to/conda/envs
mamba activate trinotate_env

N_CHUNKS=$1
MERGE_ONLY=${2:-}

OUTDIR="/path/to/project/assembly/trinotate_output"
INPUT_FILE="${OUTDIR}/Trinity.fasta.transdecoder.pep"
SPLIT_DIR="${OUTDIR}/Trinity.fasta.transdecoder.pep.split"
PARALLEL_DIR="${OUTDIR}/signalp_output_parallel"
PBS_SCRIPT="${OUTDIR}/trinotate_05_signalp6_array.pbs"

cd "${OUTDIR}"

[[ -s "${INPUT_FILE}" ]] || { echo "[FATAL] Missing input: ${INPUT_FILE}" >&2; exit 1; }
