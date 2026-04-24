#!/bin/bash

#PBS -l walltime=150:00:0
#PBS -q default@pbs-m1.metacentrum.cz
#PBS -l select=1:ncpus=30:mem=45gb:scratch_local=200gb
#PBS -N Festuca_TransDecoder_Predict
#PBS -j oe

# Author: Iris Sammarco
# Date: 03/2026
# Aim: Run TransDecoder.Predict to identify likely coding ORFs, using Diamond BLASTP (from step 01) and Pfam hmmscan (from step 02) hits as supporting evidence.
# Run: qsub trinotate_03_transdecoder.bash
# Input: Trinity.fasta, blastp.sprot.outfmt6 (step 01), pfam.domtblout (step 02)
# Output: Trinity.fasta.transdecoder.pep
#         Trinity.fasta.transdecoder.bed
#         Trinity.fasta.transdecoder.gff3
#         Trinity.fasta.transdecoder.cdna
# EggNOG must be run after this script completes, since Trinity.fasta.transdecoder.pep is generated here

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
EVAL="1e-5"

cd "${OUTDIR}"

# Per script logs:
exec 1> >(tee -a "${LOGDIR}/${PBS_JOBNAME}.log")
exec 2> >(tee -a "${LOGDIR}/${PBS_JOBNAME}.err")
echo "[INFO] Job ${PBS_JOBID} started in ${PWD}

## Sanity checks (fail-fast)
[[ -s "${ASSEMBLY}" ]] || { echo "[FATAL] Missing assembly: ${ASSEMBLY}"; exit 1; }
ln -sf "${ASSEMBLY}" Trinity.fasta

## TransDecoder.Predict - Uses BLASTP + Pfam domtblout → BETTER ORFs (20-30% fewer false positives). Predicts which ORFs are likely to be coding
if [[ ! -s "Trinity.fasta.transdecoder.pep" || ! -s ".transdecoder_predict.done" ]]; then
    echo "[INFO] TransDecoder.Predict..."
    rm -f Trinity.fasta.transdecoder.pep .transdecoder_predict.done
    TransDecoder.Predict -t Trinity.fasta \
        --retain_blastp_hits blastp.sprot.outfmt6 \
        --retain_pfam_hits pfam.domtblout \
        > "${LOGDIR}/transdecoder.predict.log" 2>&1
    touch .transdecoder_predict.done
fi
PEP_FINAL="Trinity.fasta.transdecoder.pep"

