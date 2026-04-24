#!/bin/bash

# Author: Iris Sammarco
# Date: 03/2026
# Aim: Master (non-qsub) script to split the TransDecoder predicted proteome into chunks and submit PBS array jobs for parallel signal peptide prediction with SignalP v6 (fast mode).
# Includes a merge mode to concatenate and validate chunk outputs, with duplicate ID check.
# Calls: trinotate_06_signalp6_array.pbs
# Run: bash trinotate_06_parallel_signalp6.bash 72    # split + submit
#      bash trinotate_06_parallel_signalp6.bash 72 m  # merge after all jobs finish
# Input: Trinity.fasta.transdecoder.pep (step 03)
# Output: signalp_output_parallel/signalp.merged.txt

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
PBS_SCRIPT="${OUTDIR}/trinotate_06_signalp6_array.pbs"

cd "${OUTDIR}"

[[ -s "${INPUT_FILE}" ]] || { echo "[FATAL] Missing input: ${INPUT_FILE}" >&2; exit 1; }

if [[ -z "${MERGE_ONLY}" ]]; then
    echo "[INFO] === SPLIT + SUBMIT MODE: SIGNALP (${N_CHUNKS} nominal chunks) ==="

    if [[ ! -d "${SPLIT_DIR}" || ! -f "${SPLIT_DIR}/chunk_signalp_001.pep" ]]; then
        echo "[INFO] Creating chunks in ${SPLIT_DIR}"
        rm -rf "${SPLIT_DIR}"
        mkdir -p "${SPLIT_DIR}"
        seqkit split2 -s 10000 "${INPUT_FILE}" -o "${SPLIT_DIR}" --by-size-prefix "chunk_signalp_"
    else
        echo "[INFO] Reusing existing chunks in ${SPLIT_DIR}"
    fi

    mkdir -p "${PARALLEL_DIR}" logs

    CHUNKS=$(ls "${SPLIT_DIR}"/chunk_signalp_*.pep 2>/dev/null | wc -l || echo 0)
    echo "[INFO] Found ${CHUNKS} chunk files"
    [[ ${CHUNKS} -ge 1 ]] || { echo "[FATAL] No chunks found" >&2; exit 1; }

    qsub -J 1-${CHUNKS}:1 \
        -N "signalp_parallel_${CHUNKS}" \
        -o "logs/signalp_%A_%a.out" \
        -e "logs/signalp_%A_%a.err" \
        "${PBS_SCRIPT}"

    echo "[INFO] Submitted ${CHUNKS} SignalP jobs"
    echo "[INFO] Monitor with: qstat -u \$USER | grep signalp"

elif [[ "${MERGE_ONLY}" == "m" ]]; then
    echo "[INFO] === MERGE MODE: SIGNALP ==="

    mkdir -p "${PARALLEL_DIR}" logs

    TOTAL_CHUNKS=$(ls "${SPLIT_DIR}"/chunk_signalp_*.pep 2>/dev/null | wc -l || echo 0)
    DONE_CHUNKS=$(ls "${PARALLEL_DIR}"/chunk_signalp_*.prediction_results.txt.done 2>/dev/null | wc -l || echo 0)

    echo "[INFO] Chunk completion: ${DONE_CHUNKS}/${TOTAL_CHUNKS}"
    [[ ${TOTAL_CHUNKS} -ge 1 ]] || { echo "[FATAL] No chunk inputs found in ${SPLIT_DIR}" >&2; exit 1; }
    [[ ${DONE_CHUNKS} -eq ${TOTAL_CHUNKS} ]] || { echo "[FATAL] Not all SignalP chunks are finished yet" >&2; exit 1; }

    touch .signalp_run.done
    touch .signalp.chunks.done

	echo "[INFO] Processing ${DONE_CHUNKS}/${TOTAL_CHUNKS} completed chunks..."

	SIGNALP_MERGED="${PARALLEL_DIR}/signalp.merged.txt"
    rm -f "${SIGNALP_MERGED}"
	
	# Keep the header and comments only from the first file, append only non-# lines from the others.
    first_file=1
    while read -r f; do
		[[ -s "$f" ]] || { echo "[FATAL] Missing or empty: $f" >&2; exit 1; } # fail if any file is empty

		if [[ ${first_file} -eq 1 ]]; then
			cat "$f" >> "${SIGNALP_MERGED}"
			first_file=0
		else
			grep -v '^#' "$f" >> "${SIGNALP_MERGED}"
		fi
	done < <(find "${PARALLEL_DIR}" -name "chunk_signalp_*.prediction_results.txt" | sort -V) # Guarantee correct file ordering

    [[ -s "${SIGNALP_MERGED}" ]] || { echo "[FATAL] Merged SignalP output is empty" >&2; exit 1; }

    MERGED_LINES=$(grep -vc '^#' "${SIGNALP_MERGED}" || true)
    echo "[INFO] Merged SignalP entries: ${MERGED_LINES}"

    [[ ${MERGED_LINES} -ge 1 ]] || { echo "[FATAL] No SignalP prediction lines found after merge" >&2; exit 1; }
	
	# Check for duplicate IDs - they may be created if some chunks overlap by mistake (should be empty)
    DUPLICATES=$(cut -f1 "${SIGNALP_MERGED}" | grep -v '^#' | sort | uniq -d || true)
    if [[ -n "${DUPLICATES}" ]]; then
        echo "[FATAL] Duplicate sequence IDs found in merged SignalP output:" >&2
        echo "${DUPLICATES}" >&2
        exit 1
    fi
	
	touch .signalp_convert.done .signalp.merged.done .signalp.done # Success flags

    echo "[INFO] Merged SignalP file ready: ${SIGNALP_MERGED}"
    echo "[INFO] Non-header prediction lines: ${MERGED_LINES}"
fi
