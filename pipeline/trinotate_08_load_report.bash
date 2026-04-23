#!/bin/bash

#PBS -l walltime=500:00:0
#PBS -q default@pbs-m1.metacentrum.cz
#PBS -l select=1:ncpus=2:mem=10gb:scratch_local=50gb
#PBS -N trinotate_load
#PBS -j oe

# Author: Iris Sammarco
# Date: 06/03/2026
# Aim: Load all annotation results into Trinotate SQLite DB and export the final report.
# Run: qsub trinotate_08_load_report.bash

set -euo pipefail
trap 'echo "[ERROR] Line $LINENO: ${BASH_COMMAND} failed on $(date)" >&2; exit 1' ERR

export TMPDIR=$SCRATCHDIR
module add mambaforge
export CONDA_ENVS_PATH=/path/to/conda/envs
mamba activate trinotate_env

export TRINITY_HOME="/path/to/conda/envs/trinotate_env/opt/trinity-2.8.5"
export PATH="/path/to/project/tmhmm-2.0c/bin:$PATH"

OUTDIR="/path/to/project/assembly/trinotate_output"
LOGDIR="${OUTDIR}/logs"
PREFIX="Festuca_rubra"

mkdir -p "${OUTDIR}" "${LOGDIR}"

exec > >(tee -a "${LOGDIR}/${PBS_JOBNAME}.log") 2>&1

echo "[INFO] Job ${PBS_JOBID} started"
echo "[INFO] Using SCRATCHDIR: $SCRATCHDIR"
