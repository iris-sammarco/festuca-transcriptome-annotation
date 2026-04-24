#!/bin/bash

#PBS -l walltime=150:00:0
#PBS -q default@pbs-m1.metacentrum.cz
#PBS -l select=1:ncpus=12:mem=45gb:scratch_local=200gb
#PBS -N Festuca_init_tmhmm
#PBS -j oe

# Author: Iris Sammarco
# Date: 03/2026
# Aim: Predict transmembrane helices in the TransDecoder-predicted proteome using TMHMM v2.
# Requires Trinity.fasta.transdecoder.pep from step 03.
# Run: qsub trinotate_05_tmhmm.bash
# Input: Trinity.fasta.transdecoder.pep (step 03)
# Output: tmhmm.out

set -euo pipefail
trap 'echo "[ERROR] Line $LINENO: ${BASH_COMMAND} failed on $(date)" >&2; exit 1' ERR

export TMPDIR=$SCRATCHDIR
module add mambaforge
export CONDA_ENVS_PATH=/path/to/conda/envs
mamba activate trinotate_env
export TRINITY_HOME="/path/to/conda/envs/trinotate_env/opt/trinity-2.8.5"
export PATH="/path/to/project/tmhmm-2.0c/bin:$PATH"
export TRINOTATE_DATA_DIR="/path/to/project/assembly/trinotate_data"
export EGGNOG_DATA_DIR="${TRINOTATE_DATA_DIR}/eggnog_data"

ASSEMBLY="/path/to/project/assembly/Trinity_output.Trinity.fasta"
OUTDIR="/path/to/project/assembly/trinotate_output"
LOGDIR="${OUTDIR}/logs"

mkdir -p "${OUTDIR}" "${LOGDIR}"
cd "${OUTDIR}"

# Per script logs:
exec 1> >(tee -a "${LOGDIR}/${PBS_JOBNAME}.log")
exec 2> >(tee -a "${LOGDIR}/${PBS_JOBNAME}.err")
echo "[INFO] Job ${PBS_JOBID} started in ${PWD}"

## Sanity checks (fail-fast)
[[ -s "${ASSEMBLY}" ]] || { echo "[FATAL] Missing assembly: ${ASSEMBLY}"; exit 1; }
ln -sf "${ASSEMBLY}" Trinity.fasta

## TMHMM - predict transmembrane domains
PEP_FINAL="Trinity.fasta.transdecoder.pep"
if [[ ! -s ".tmhmm.done" ]]; then
    echo "[INFO] TMHMM predictions..."
    rm -f tmhmm.out .tmhmm.done
    
    tmhmm --short < "${PEP_FINAL}" > tmhmm.out 2>"${LOGDIR}/tmhmm.log" \
        || { echo "[FATAL] TMHMM failed. Verify PATH."; exit 1; }

    touch .tmhmm.done
fi
