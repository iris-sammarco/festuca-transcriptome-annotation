#!/bin/bash

#PBS -l walltime=150:00:0
#PBS -q default@pbs-m1.metacentrum.cz
#PBS -l select=1:ncpus=28:mem=40gb:scratch_local=200gb
#PBS -N Festuca_Main
#PBS -j oe

# Author: Iris Sammarco
# Date: 06/03/2026
# Aim: Generate the Trinity gene-to-transcript map, build the Diamond SwissProt database (if absent), extract candidate long ORFs with TransDecoder.LongOrfs, and run Diamond BLASTP against SwissProt to provide homology evidence for downstream ORF prediction (step 03).
# Run: qsub trinotate_01_longorfs_blastp.bash
# Input: Trinity_output.Trinity.fasta, uniprot_sprot.pep
# Output: Trinity.fasta.gene_trans_map
#         Trinity.fasta.transdecoder_dir/longest_orfs.pep (input for Pfam, step 02)
#         blastp.sprot.outfmt6

set -euo pipefail  # Exit on error, undefined vars, pipe failures
trap 'echo "[ERROR] Line $LINENO: ${BASH_COMMAND} failed on $(date)" >&2; exit 1' ERR # debug

# Environment & Paths
export TMPDIR=$SCRATCHDIR # This will force the application to place the temporary files into scratch directory instead of a /tmp directory
module add mambaforge
export CONDA_ENVS_PATH=/path/to/conda/envs
mamba activate trinotate_env
export TRINITY_HOME="/path/to/conda/envs/trinotate_env/opt/trinity-2.8.5"
export PATH="/path/to/project/tmhmm-2.0c/bin:$PATH"
export TRINOTATE_DATA_DIR="/path/to/project/assembly/trinotate_data" # expects the files: Trinotate.sqlite and uniprot_sprot.fa
export EGGNOG_DATA_DIR="${TRINOTATE_DATA_DIR}/eggnog_data"

# Config
ASSEMBLY="/path/to/project/assembly/Trinity_output.Trinity.fasta"
PREFIX="Festuca_rubra"
OUTDIR="/path/to/project/assembly/trinotate_output"
LOGDIR="${OUTDIR}/logs"
THREADS=28 # Must match ncpus in PBS header
EVAL="1e-5"  # Strict Diamond E-value
MAX_TARGETS=5  # Maximum number of target sequences to report alignments for (default=25)

mkdir -p "${OUTDIR}" "${LOGDIR}"
cd "${OUTDIR}"

# Per script logs:
exec 1> >(tee -a "${LOGDIR}/${PBS_JOBNAME}.log")
exec 2> >(tee -a "${LOGDIR}/${PBS_JOBNAME}.err")
echo "[INFO] Job ${PBS_JOBID} started in ${PWD} | Threads: ${THREADS}"

# Create symlinks of the db files in the OUTDIR:
ln -sf "${TRINOTATE_DATA_DIR}"/Trinotate.* .
ln -sf "${TRINOTATE_DATA_DIR}"/Pfam-A.hmm* .
ln -sf "${TRINOTATE_DATA_DIR}"/uniprot_sprot.pep .
ln -sf "${TRINOTATE_DATA_DIR}"/{uniprot_sprot.diamond,Rfam.cm,Rfam.clanin} . 2>/dev/null || true

## Sanity checks (fail-fast)
[[ -s "${ASSEMBLY}" ]] || { echo "[FATAL] Missing assembly: ${ASSEMBLY}"; exit 1; }
ln -sf "${ASSEMBLY}" Trinity.fasta

## Step 1: Build Diamond DB if missing
DB="${TRINOTATE_DATA_DIR}/uniprot_sprot.diamond.dmnd"
if [[ ! -s "${DB}" || ! -s ".diamond_db.done" ]]; then # Checks if output file exists AND has size >0 (-s). AND .step.done missing -> Runs commands.
    echo "[INFO] Building Diamond SwissProt DB..."
    diamond makedb --in "${TRINOTATE_DATA_DIR}/uniprot_sprot.pep" -d "${TRINOTATE_DATA_DIR}/uniprot_sprot.diamond" --threads ${THREADS} > "${LOGDIR}/diamond_makedb.log" 2>&1
    touch .diamond_db.done
fi
[[ -s "${DB}" ]] || { echo "[FATAL] Diamond DB build failed"; exit 1; }

## Step 2: Trinity gene→transcript map
if [[ ! -s "Trinity.fasta.gene_trans_map" || ! -s ".gene_trans_map.done" ]]; then
    echo "[INFO] Generating gene-trans map..."
    rm -f Trinity.fasta.gene_trans_map .gene_trans_map.done
    ${TRINITY_HOME}/util/support_scripts/get_Trinity_gene_to_trans_map.pl Trinity.fasta > Trinity.fasta.gene_trans_map 2> "${LOGDIR}/gene_trans_map.log"
    touch .gene_trans_map.done
fi
[[ -s "Trinity.fasta.gene_trans_map" ]] || { echo "[FATAL] Gene-trans map failed"; exit 1; }

## Step 3: TransDecoder LongOrfs (candidate long ORFs)
if [[ ! -s "Trinity.fasta.transdecoder_dir/longest_orfs.pep" || ! -s ".transdecoder_longorfs.done" ]]; then
    echo "[INFO] TransDecoder.LongOrfs..."
    rm -rf Trinity.fasta.transdecoder_dir .transdecoder_longorfs.done
    TransDecoder.LongOrfs -t Trinity.fasta > "${LOGDIR}/transdecoder.longorfs.log" 2>&1
    touch .transdecoder_longorfs.done
fi
PEP_LONG="Trinity.fasta.transdecoder_dir/longest_orfs.pep"

## Step 4: Direct Diamond SwissProt BLASTP (for TransDecoder evidence)
if [[ ! -s "blastp.sprot.outfmt6" || ! -s ".blastp_evidence.done" ]]; then
    echo "[INFO] Diamond BLASTP (TransDecoder evidence)..."
    rm -f blastp.sprot.outfmt6 .blastp_evidence.done
    diamond blastp \
        --db "${TRINOTATE_DATA_DIR}/uniprot_sprot.diamond" \
        --query "${PEP_LONG}" \
        --out blastp.sprot.outfmt6 \
        --evalue ${EVAL} --max-target-seqs ${MAX_TARGETS} --threads ${THREADS} --outfmt 6 \
        > "${LOGDIR}/blastp.evidence.log" 2>&1 # --outfmt 6 controls the output format (BLAST tabular)
    touch .blastp_evidence.done
fi
