#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "Usage: $0 <mopac_exe> <template_mop> <out_dir> <threads_csv> [runs]" >&2
  echo "Example: $0 ./build/mopac receptor_mozyme_pm6_1scf_profile.mop /tmp/bench 1,2,4,8,16 3" >&2
  exit 1
fi

mopac_exe=$1
template_mop=$2
out_dir=$3
threads_csv=$4
runs=${5:-3}

mkdir -p "$out_dir"
IFS=',' read -r -a threads <<< "$threads_csv"

for t in "${threads[@]}"; do
  for r in $(seq 1 "$runs"); do
    base="$out_dir/receptor_diag2_autotune_t${t}_r${r}"
    mop="${base}.mop"
    cp "$template_mop" "$mop"
    sed -i -E "1s/THREADS=[0-9]+/THREADS=${t}/" "$mop"
    /usr/bin/time -f "%e" -o "${base}.time" \
      env OMP_NUM_THREADS="$t" OPENBLAS_NUM_THREADS=1 "$mopac_exe" "$mop" \
      > "${base}.stdout" 2>&1
  done
done

echo "Benchmark completed in: $out_dir"
