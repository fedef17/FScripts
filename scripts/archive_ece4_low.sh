#!/usr/bin/env bash
#SBATCH --job-name=archive_ece4
#SBATCH --output=archive_ece4_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
 
set -euo pipefail
 
# --- config ---
ORIG_DIR="${SCRATCH}/ece4/"
DEST_DIR="ec:/ccff/ece4/tuning/TL63/"

if [[ $# -lt 1 ]]; then
  echo "Usage: $(basename "$0") <name> [name2 ...]"
  exit 1
fi
NAMES=("$@")

#NAMES=(
#ll00 ll01 ll02 ll03 ll04 ll05 ll06 ll07 ll08 ll09 ll10 ll11 ll12 ll13 ll14 ll15 ll16 ll17 ll18 ll19 ll20 ll21 ll22 ll23 ll24 ll25 ll26 ll27 ll28 ll29 ll30 ll31 ll32 ll33 ll34 ll35 ll36 ll37 ll38 ll39 ll40 ll41 ll42 ll57 ll58 lltu llv4 llv5 llv7
#)
 
ARCHIVE_SUBDIRS=(
  "ecmean"
  "restart"
  "output"
)
# --------------

# Helper: create archive in orig_dir, okcp it to dest, then remove local copy
transfer() {
  local archive="$1"   # full path in orig_dir (temp location)
  local dest="$2"      # destination path (dir or file)
  echo "[$(date '+%H:%M:%S')]     ecp '$archive' -> '$dest'"
  ecp "$archive" "$dest"
  rm -f "$archive"
}
 
for name in "${NAMES[@]}"; do
  src="${ORIG_DIR}/${name}"
 
  if [[ ! -d "$src" || -L "$src" ]]; then
    echo "[SKIP] '$src' is not a real directory."
    continue
  fi
 
  echo "[$(date '+%H:%M:%S')] Processing '$name' ..."
 
  # 1. Recreate first-level structure under dest
  for subdir in "${ARCHIVE_SUBDIRS[@]}"; do
    emkdir -p "${DEST_DIR}/${name}/${subdir}"
  done
 
  # 2. Archive regular subdirs (one archive per subdir, skip 'output')
  for subdir in "${ARCHIVE_SUBDIRS[@]}"; do
    [[ "$subdir" == "output" ]] && continue
 
    subdir_src="${src}/${subdir}"
 
    if [[ ! -d "$subdir_src" || -L "$subdir_src" ]]; then
      echo "[SKIP] '$subdir_src' is not a real directory."
      continue
    fi
 
    archive_name="${name}_${subdir}.tar.gz"
    archive_tmp="${src}/${archive_name}"
    echo "[$(date '+%H:%M:%S')]   Archiving '$subdir_src' ..."
    tar -czf "$archive_tmp" -C "$subdir_src" .
    transfer "$archive_tmp" "${DEST_DIR}/${name}/${subdir}/${archive_name}"
  done
 
  # 3. Handle 'output' subdir: group files by decade
  output_src="${src}/output"
 
  if [[ ! -d "$output_src" || -L "$output_src" ]]; then
    echo "[SKIP] '$output_src' is not a real directory."
    continue
  fi

  for model in oifs nemo; do 
    # Collect files by decade from filenames like: name_*_1990-1999.nc
    declare -A decade_files
    while IFS= read -r -d '' f; do
      fname=$(basename "$f")
      if [[ "$fname" =~ _([0-9]{4})-[0-9]{4}\.nc$ ]]; then
        year="${BASH_REMATCH[1]}"
        decade=$(( (year / 10) * 10 ))
        decade_files[$decade]+="$fname"$'\n'
      fi
    done < <(find "${output_src}/${model}/" -maxdepth 1 -type f -not -type l -name "*.nc" -print0)
 
    for decade in $(echo "${!decade_files[@]}" | tr ' ' '\n' | sort -n); do
      decade_end=$(( decade + 9 ))
      archive_name="output_${model}_${decade}-${decade_end}.tar.gz"
      archive_tmp="${src}/${archive_name}"
      echo "[$(date '+%H:%M:%S')]   Archiving output $model decade ${decade}-${decade_end} ..."
 
      tmp_list=$(mktemp)
      echo -n "${decade_files[$decade]}" > "$tmp_list"
      tar -czf "$archive_tmp" -C "$output_src/${model}" --files-from="$tmp_list"
      rm -f "$tmp_list"
 
      transfer "$archive_tmp" "${DEST_DIR}/${name}/output/${archive_name}"
    done
 
    unset decade_files
    echo "[$(date '+%H:%M:%S')] Done: $name"
  done
done
