#!/bin/bash
#SBATCH --job-name=rec_tape   # Job name
#SBATCH --output=rec_tape_%j.out # Output file
#SBATCH --error=rec_tape_%j.err  # Error file
#SBATCH --time=48:00:00                # Wall time limit
#SBATCH --nodes=1                      # Number of nodes
#SBATCH --ntasks-per-node=1            # Number of tasks
#SBATCH --mem=8G                       # Memory per node

# Validate arguments
if [ $# -ne 3 ]; then
    echo "ERROR: Required arguments missing"
    echo "Usage: $0 <expname> <start_year> <end_year>"
    exit 1
fi

expname=$1
START_YEAR=$2
END_YEAR=$3

# Configuration
#BASE_DIR="ec:/ccff/ece3/tunecs/${expname}/cmorized/"      # Directory containing yearly tar archives (full path on tape)
BASE_DIR="ec:/ccpd/${expname}/cmorized/"      # Directory containing yearly tar archives (full path on tape)
ARCHIVE_PATTERN="${expname}_cmorized_%Y.part.??" # Archive naming pattern (%Y will be replaced with year)
OUTPUT_DIR="$SCRATCH/tunecs_coupled/${expname}/"   # Final output directory

#VARIABLES=("ts" "tas" "hus" "ta" "rsus" "rsds" "rlut" "rsut" "rsdt" "rlutcs" "rsutcs") # Variables to keep
#VARIABLES=("tas" "rlut" "rsut" "rsutcs" "rlutcs") # Variables to keep
#VARIABLES=("hur")
VARIABLES=("clt" "clivi" "clwvi") # Variables to keep
miptab='Amon' # to keep all tabs, write "*" here

echo "Starting job $SLURM_JOB_ID at $(date)"

# Create output directory structure
mkdir -p "$OUTPUT_DIR"
for var in "${VARIABLES[@]}"; do
    mkdir -p "$OUTPUT_DIR/$var"
done

cd $OUTPUT_DIR

# Process each year
for (( YEAR=START_YEAR; YEAR<=END_YEAR; YEAR++ )); do
    echo "[$(date)] Processing year: $YEAR"
    
    okfile=${ARCHIVE_PATTERN//%Y/$YEAR}

    # Construct paths
    YEAR_FOLDER="${BASE_DIR}/cmor_${YEAR}"
    ARCHIVE_FILE="${YEAR_FOLDER}/${okfile}"
    
    # Create temporary directory for this year
    TEMP_YEAR_DIR=$(mktemp -d -p "$OUTPUT_DIR" temp_${YEAR}_XXXXXX)
    
    # Copy all split parts matching the pattern
    echo "[$(date)] Copying archive parts for year $YEAR..."
    ecp "${ARCHIVE_FILE}" .
    if [[ $? -ne 0 ]]; then
        echo "ERROR: Failed to copy archive parts for year $YEAR"
        continue
    fi

    # Check how many parts are on disk
    PART_FILES=( $(ls ${okfile} 2>/dev/null | sort) )
    N_PARTS=${#PART_FILES[@]}
    if [[ $N_PARTS -eq 0 ]]; then
        echo "ERROR: No archive parts found on disk for year $YEAR after copy"
        continue
    fi
    echo "[$(date)] Found $N_PARTS part(s) on disk: ${PART_FILES[*]}"

    # Rebuild the split archive if needed, or use single part directly
    if [[ $N_PARTS -gt 1 ]]; then
        RECONSTRUCTED="cmorized_${YEAR}_reconstructed.tar"
        echo "[$(date)] Reconstructing archive from $N_PARTS parts..."
        cat "${PART_FILES[@]}" > "$RECONSTRUCTED"
        if [[ $? -ne 0 ]]; then
            echo "ERROR: Failed to reconstruct archive for year $YEAR"
            continue
        fi
    else
        RECONSTRUCTED="${PART_FILES[0]}"
        echo "[$(date)] Single part found"
    fi

    tar -xf "$RECONSTRUCTED" -C "$TEMP_YEAR_DIR"
    
    # Process extracted files
    cd "$TEMP_YEAR_DIR" || { echo "Failed to cd to $TEMP_YEAR_DIR" >&2; exit 1; }
    
    # Find and copy specified variables
    for var in "${VARIABLES[@]}"; do
        echo "[$(date)] Processing variable $var for year $YEAR"
        find "${TEMP_YEAR_DIR}" -type f -name "*${var}_${miptab}*" -exec mv {} "${OUTPUT_DIR}/${var}/" \;
    done
    
    # Clean up
    cd - || exit
    rm -rf "$TEMP_YEAR_DIR"
    rm -f "${PART_FILES[@]}" "$RECONSTRUCTED"
done

echo "[$(date)] Processing complete. Results are in $OUTPUT_DIR"
echo "Job $SLURM_JOB_ID finished at $(date)"
