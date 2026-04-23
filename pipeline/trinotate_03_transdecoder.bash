#!/bin/bash

#PBS -l walltime=150:00:0
#PBS -q default@pbs-m1.metacentrum.cz
#PBS -l select=1:ncpus=30:mem=45gb:scratch_local=200gb
#PBS -N Festuca_TransDecoder_Predict
#PBS -j oe

# Author: Iris Sammarco
# Date: 06/03/2026
# Aim: TransDecoder ORF prediction on Trinity assembly.
# Run: qsub trinotate_03_transdecoder.bash

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
PREFIX="Festuca_rubra"
OUTDIR="/path/to/project/assembly/trinotate_output"
LOGDIR="${OUTDIR}/logs"
THREADS=28
EVAL="1e-5"
MAX_TARGETS=5

cd "${OUTDIR}"

# Per script logs:
exec 1> >(tee -a "${LOGDIR}/${PBS_JOBNAME}.log")
exec 2> >(tee -a "${LOGDIR}/${PBS_JOBNAME}.err")
echo "[INFO] Job ${PBS_JOBID} started in ${PWD} | Threads: ${THREADS}"

## Sanity checks (fail-fast)
[[ -s "${ASSEMBLY}" ]] || { echo "[FATAL] Missing assembly: ${ASSEMBLY}"; exit 1; }
ln -sf "${ASSEMBLY}" Trinity.fasta

## Step 5: TransDecoder.Predict - Uses BLASTP + Pfam domtblout → BETTER ORFs (20-30% fewer false positives). Predicts which ORFs are likely to be coding
if [[ ! -s "Trinity.fasta.transdecoder.pep" || ! -s ".transdecoder_predict.done" ]]; then
    echo "[INFO] TransDecoder.Predict..."
    rm -f Trinity.fasta.transdecoder.pep .transdecoder_predict.done
    TransDecoder.Predict -t Trinity.fasta \
        --retain_blastp_hits blastp.sprot.outfmt6 \ # produced by step 01
        --retain_pfam_hits pfam.domtblout \ # produced by step 02
        > "${LOGDIR}/transdecoder.predict.log" 2>&1
    touch .transdecoder_predict.done
fi
PEP_FINAL="Trinity.fasta.transdecoder.pep"

