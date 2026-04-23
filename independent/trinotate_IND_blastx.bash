#!/bin/bash

#PBS -l walltime=48:00:0
#PBS -q default@pbs-m1.metacentrum.cz
#PBS -l select=1:ncpus=28:mem=40gb:scratch_local=200gb
#PBS -N Festuca_BlastX
#PBS -j oe

# Author: Iris Sammarco
# Date: 06/03/2026
# Aim: Run the BLASTX step independently from the main pipeline.
#      This step is not position-dependent and can be run at any time before the final load (step 08).
# Run: qsub trinotate_IND_blastx.bash

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

mkdir -p "${OUTDIR}" "${LOGDIR}"
cd "${OUTDIR}"

# Per script logs:
exec 1> >(tee -a "${LOGDIR}/${PBS_JOBNAME}.log")
exec 2> >(tee -a "${LOGDIR}/${PBS_JOBNAME}.err")
echo "[INFO] Job ${PBS_JOBID} started in ${PWD} | Threads: ${THREADS}"

# Create symlinks of the db files in the OUTDIR:
ln -sf "${TRINOTATE_DATA_DIR}"/Trinotate.sqlite* .
ln -sf "${TRINOTATE_DATA_DIR}"/Pfam-A.hmm* .
ln -sf "${TRINOTATE_DATA_DIR}"/uniprot_sprot.pep .
ln -sf "${TRINOTATE_DATA_DIR}"/{uniprot_sprot.diamond,Rfam.cm,Rfam.clanin} . 2>/dev/null || true

## Sanity checks (fail-fast)
[[ -s "${ASSEMBLY}" ]] || { echo "[FATAL] Missing assembly: ${ASSEMBLY}"; exit 1; }
ln -sf "${ASSEMBLY}" Trinity.fasta

echo "[INFO] BLASTX-only job started $(date)"

## STEP 7: Diamond BLASTX (full transcripts vs SwissProt)
if [[ ! -s "blastx.sprot.outfmt6" || ! -s ".blastx.done" ]]; then
    echo "[INFO] Diamond BLASTX..."
    rm -f blastx.sprot.outfmt6 .blastx.done
    diamond blastx \
        --db "${TRINOTATE_DATA_DIR}/uniprot_sprot.diamond" \
        --query Trinity.fasta \
        --out blastx.sprot.outfmt6 \
        --evalue ${EVAL} --max-target-seqs ${MAX_TARGETS} --threads ${THREADS} --outfmt 6 \
        > "${LOGDIR}/blastx.log" 2>&1
    touch .blastx.done
	echo "[SUCCESS] BLASTX complete $(date)"
fi
