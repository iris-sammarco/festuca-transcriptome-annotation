#!/bin/bash

#PBS -l walltime=150:00:0
#PBS -q default@pbs-m1.metacentrum.cz
#PBS -l select=1:ncpus=12:mem=45gb:scratch_local=200gb
#PBS -N Festuca_init_tmhmm
#PBS -j oe

# Author: Iris Sammarco
# Date: 06/03/2026
# Aim: Annotate Festuca rubra Trinity assembly with Trinotate + rice/wheat curation as no reference genome is available for F. rubra. When each step terminates successfully, it creates an output file with the "done" flag, so if the script needs to be resumed it will skip all the steps that have the files complete.
# Separate the longest steps in separate scripts.
# Outputs: annotation_report.xls with Pfam domains, TM/signal peptides, homology
# Run: qsub Trinotate_main.bash

set -euo pipefail  # Exit on error, undefined vars, pipe failures
trap 'echo "[ERROR] Line $LINENO: ${BASH_COMMAND} failed on $(date)" >&2; exit 1' ERR # debug

# Environment & Paths
export TMPDIR=$SCRATCHDIR # This will force the application to place the temporary files into scratch directory instead of a /tmp directory
module add mambaforge
export CONDA_ENVS_PATH=/storage/pruhonice1-ibot/home/irissammarco/.conda/envs
mamba activate trinotate_env
export TRINITY_HOME="/storage/pruhonice1-ibot/home/irissammarco/.conda/envs/trinotate_env/opt/trinity-2.8.5"
export PATH="/storage/plzen1/home/irissammarco/Festuca_RNA_assembly/tmhmm-2.0c/bin:$PATH"
export TRINOTATE_DATA_DIR="/storage/plzen1/home/irissammarco/Festuca_RNA_assembly/assembly/trinotate_data" # expects the files: Trinotate.sqlite and uniprot_sprot.fa
export EGGNOG_DATA_DIR="${TRINOTATE_DATA_DIR}/eggnog_data"

# Config
ASSEMBLY="/storage/plzen1/home/irissammarco/Festuca_RNA_assembly/assembly/Trinity_output.Trinity.fasta"
PREFIX="Festuca_rubra"
OUTDIR="/storage/plzen1/home/irissammarco/Festuca_RNA_assembly/assembly/trinotate_output"
LOGDIR="${OUTDIR}/logs"
THREADS=12
EVAL="1e-5"  # Strict Diamond E-value
MAX_TARGETS=5  # Maximum number of target sequences to report alignments for (default=25)

mkdir -p "${OUTDIR}" "${LOGDIR}"
cd "${OUTDIR}"

# Per script logs:
exec 1> >(tee -a "${LOGDIR}/${PBS_JOBNAME}.log")
exec 2> >(tee -a "${LOGDIR}/${PBS_JOBNAME}.err")
echo "[INFO] Job ${PBS_JOBID} started in ${PWD} | Threads: ${THREADS}"

## Sanity checks (fail-fast)
[[ -s "${ASSEMBLY}" ]] || { echo "[FATAL] Missing assembly: ${ASSEMBLY}"; exit 1; }
ln -sf "${ASSEMBLY}" Trinity.fasta

## STEP 9: TMHMM - predict transmembrane domains
PEP_FINAL="Trinity.fasta.transdecoder.pep"
if [[ ! -s ".tmhmm.done" ]]; then
    echo "[INFO] TMHMM predictions..."
    rm -f tmhmm.out .tmhmm.done
    
    tmhmm --short < "${PEP_FINAL}" > tmhmm.out 2>"${LOGDIR}/tmhmm.log" \
        || { echo "[FATAL] TMHMM failed. Verify PATH."; exit 1; }

    touch .tmhmm.done
fi
