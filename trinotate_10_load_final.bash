#!/bin/bash

#PBS -l walltime=500:00:0
#PBS -q default@pbs-m1.metacentrum.cz
#PBS -l select=1:ncpus=2:mem=10gb:scratch_local=50gb
#PBS -N trinotate_load
#PBS -j oe

# Author: Iris Sammarco
# Date: 06/03/2026
# Aim: Run the final step to load all data together.
# Run: qsub trinotate_10_load_final.bash

set -euo pipefail  # Exit on error, undefined vars, pipe failures
trap 'echo "[ERROR] Line $LINENO: ${BASH_COMMAND} failed on $(date)" >&2; exit 1' ERR # debug

# Environment & Paths
export TMPDIR=$SCRATCHDIR # This will force the application to place the temporary files into scratch directory instead of a /tmp directory
module add mambaforge
export CONDA_ENVS_PATH=/storage/pruhonice1-ibot/home/irissammarco/.conda/envs
mamba activate trinotate_env

export TRINITY_HOME="/storage/pruhonice1-ibot/home/irissammarco/.conda/envs/trinotate_env/opt/trinity-2.8.5"
export PATH="/storage/plzen1/home/irissammarco/Festuca_RNA_assembly/tmhmm-2.0c/bin:$PATH"

# Config
OUTDIR="/storage/plzen1/home/irissammarco/Festuca_RNA_assembly/assembly/trinotate_output"
LOGDIR="${OUTDIR}/logs"
PREFIX="Festuca_rubra"

mkdir -p "${OUTDIR}" "${LOGDIR}"

# Logging
exec > >(tee -a "${LOGDIR}/${PBS_JOBNAME}.log") 2>&1

echo "[INFO] Job ${PBS_JOBID} started"
echo "[INFO] Using SCRATCHDIR: $SCRATCHDIR"

# Copy everything to SCRATCH, critical for accelerating the Job
echo "[INFO] Copying data to scratch..."

cp ${OUTDIR}/Trinity.fasta $SCRATCHDIR/
cp ${OUTDIR}/Trinity.fasta.gene_trans_map $SCRATCHDIR/
cp ${OUTDIR}/Trinity.fasta.transdecoder.pep $SCRATCHDIR/
cp ${OUTDIR}/Trinotate.sqlite $SCRATCHDIR/
cp ${OUTDIR}/blastp.rice.outfmt6 $SCRATCHDIR/
cp ${OUTDIR}/blastp.wheat.outfmt6 $SCRATCHDIR/
cp ${OUTDIR}/blastx.sprot.outfmt6 $SCRATCHDIR/
cp ${OUTDIR}/blastp.sprot.outfmt6 $SCRATCHDIR/
cp ${OUTDIR}/pfam.domtblout $SCRATCHDIR/
cp ${OUTDIR}/tmhmm.out $SCRATCHDIR/
cp ${OUTDIR}/rfam.tblout $SCRATCHDIR/ || true

cp -r ${OUTDIR}/signalp_output_parallel $SCRATCHDIR/

cd $SCRATCHDIR

# I'll init a clean database to avoid duplicate rows:
echo "[INFO] Initializing clean Trinotate DB..."

Trinotate Trinotate.sqlite init \
  --gene_trans_map Trinity.fasta.gene_trans_map \
  --transcript_fasta Trinity.fasta \
  --transdecoder_pep Trinity.fasta.transdecoder.pep

# Safety check: Fail if missing critical files
for f in blastx.sprot.outfmt6 pfam.domtblout tmhmm.out signalp_output_parallel/signalp.merged.txt; do
    [[ -s "$f" ]] || { echo "[FATAL] Missing $f"; exit 1; }
done

## Load Step
rm -f .trinotate.loads.done

echo "[INFO] LOAD step started $(date)"

Trinotate Trinotate.sqlite LOAD_custom_blast --outfmt6 blastp.rice.outfmt6 --prog blastp --dbtype rice_proteome
Trinotate Trinotate.sqlite LOAD_custom_blast --outfmt6 blastp.wheat.outfmt6 --prog blastp --dbtype wheat_proteome
Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.sprot.outfmt6
Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.sprot.outfmt6
Trinotate Trinotate.sqlite LOAD_pfam pfam.domtblout
Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out
Trinotate Trinotate.sqlite LOAD_signalp signalp_output_parallel/signalp.merged.txt

#Trinotate Trinotate.sqlite LOAD_eggnog eggnog.emapper.annotations 2> "${LOGDIR}/load.eggnog.log" 2>&1 # the version of Trinotate I'm using (3.2.2) doesnt support this load
#Trinotate Trinotate.sqlite LOAD_rfam rfam.tblout 2> "${LOGDIR}/load.rfam.log" 2>&1 # the version of Trinotate I'm using (3.2.2) doesnt support this load
touch .trinotate.loads.done

echo "[INFO] LOAD completed $(date)"

## Report step
#rm -f ${PREFIX}.annotation_report*.xls
echo "[INFO] REPORT step started $(date)"

Trinotate Trinotate.sqlite report -E 1e-20 \
	--pfam_cutoff DGC | gzip > "${PREFIX}.annotation_report_strict.xls.gz"

echo "[INFO] REPORT completed $(date)"

## Copy results back
echo "[INFO] Copying results back..."

cp ${PREFIX}.annotation_report_strict.xls.gz ${OUTDIR}/

echo "✅ DONE"