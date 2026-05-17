#!/bin/bash
# Parallel scan runner for old program exploration
# Runs gamma, bound, weight, operator, and fine compression scans
set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate analyse

PROJDIR="/mnt/d/Codes/airesearch/NewINCGengrate"
OUTDIR="$PROJDIR/old_program_exploration"
RAWDIR="$OUTDIR/raw"
LOGDIR="$OUTDIR/logs"
cd "$PROJDIR"

MAX_PARALLEL=6  # 13650HX has 14 cores, use 6 parallel runs

# Track running jobs
declare -a PIDS=()
declare -a NAMES=()

wait_for_slot() {
  while [ ${#PIDS[@]} -ge $MAX_PARALLEL ]; do
    NEW_PIDS=()
    NEW_NAMES=()
    for i in "${!PIDS[@]}"; do
      if kill -0 "${PIDS[$i]}" 2>/dev/null; then
        NEW_PIDS+=("${PIDS[$i]}")
        NEW_NAMES+=("${NAMES[$i]}")
      else
        wait "${PIDS[$i]}" 2>/dev/null
        echo "[$(date +%H:%M:%S)] DONE: ${NAMES[$i]}"
      fi
    done
    PIDS=("${NEW_PIDS[@]+"${NEW_PIDS[@]}"}")
    NAMES=("${NEW_NAMES[@]+"${NEW_NAMES[@]}"}")
    if [ ${#PIDS[@]} -ge $MAX_PARALLEL ]; then
      sleep 5
    fi
  done
}

launch_run() {
  local name="$1"
  shift
  wait_for_slot
  echo "[$(date +%H:%M:%S)] START: $name"
  "$@" > "$LOGDIR/${name}.log" 2>&1 &
  PIDS+=($!)
  NAMES+=("$name")
}

wait_all() {
  for i in "${!PIDS[@]}"; do
    wait "${PIDS[$i]}" 2>/dev/null
    echo "[$(date +%H:%M:%S)] DONE: ${NAMES[$i]}"
  done
  PIDS=()
  NAMES=()
}

# ===========================================================
# SECTION C: v18 Gamma Scan (skip gamma=0, already done)
# ===========================================================
echo "=== SECTION C: v18 Gamma Scan ==="
for gamma in 0.02 0.05 0.10 0.25 0.50 1.0; do
  OUT_SUB="$RAWDIR/v18_poly_gamma_${gamma}"
  mkdir -p "$OUT_SUB"
  launch_run "v18_poly_gamma_${gamma}" bash -c "
    cd '$OUT_SUB'
    '$PROJDIR/inc_ricci_v18' --nx 16 --ny 16 \
      --radius_seed 101 --target_dist lognormal --target_delta 0.25 \
      --target_radius_ratio 5 --coord_mode soft_relax --target_phi 0.84 \
      --relax_steps 20000 --use_voropp 1 --compute_nm_theta 1 \
      --decompose_theta 1 --apply_to_polydisperse 1 \
      --cycles 200 --shape_steps_per_cycle 50 --ricci_steps_per_cycle 50 \
      --shape_gamma ${gamma} --bounded_du 0.10 --progress_every 50
    echo EXIT=\$?
  "
done

# ===========================================================
# SECTION D: v18 Bounded_du Scan
# ===========================================================
echo "=== SECTION D: v18 Bounded_du Scan ==="
for bdu in 0.03 0.05 0.15 0.25 0.40; do
  # 0.10 already done in gamma scan baseline
  OUT_SUB="$RAWDIR/v18_poly_bound_${bdu}"
  mkdir -p "$OUT_SUB"
  launch_run "v18_poly_bound_${bdu}" bash -c "
    cd '$OUT_SUB'
    '$PROJDIR/inc_ricci_v18' --nx 16 --ny 16 \
      --radius_seed 101 --target_dist lognormal --target_delta 0.25 \
      --target_radius_ratio 5 --coord_mode soft_relax --target_phi 0.84 \
      --relax_steps 20000 --use_voropp 1 --compute_nm_theta 1 \
      --decompose_theta 1 --apply_to_polydisperse 1 \
      --cycles 200 --shape_steps_per_cycle 20 --ricci_steps_per_cycle 100 \
      --shape_gamma 0.05 --bounded_du ${bdu} --progress_every 50
    echo EXIT=\$?
  "
done

# ===========================================================
# SECTION E: v18 Weight Scan
# ===========================================================
echo "=== SECTION E: v18 Weight Scan ==="

# shape_biased
OUT_SUB="$RAWDIR/v18_weight_shape_biased"
mkdir -p "$OUT_SUB"
launch_run "v18_weight_shape_biased" bash -c "
  cd '$OUT_SUB'
  '$PROJDIR/inc_ricci_v18' --nx 16 --ny 16 \
    --radius_seed 101 --target_dist lognormal --target_delta 0.25 \
    --target_radius_ratio 5 --coord_mode soft_relax --target_phi 0.84 \
    --relax_steps 20000 --use_voropp 1 --compute_nm_theta 1 \
    --decompose_theta 1 --apply_to_polydisperse 1 \
    --cycles 200 --shape_steps_per_cycle 50 --ricci_steps_per_cycle 50 \
    --shape_gamma 0.05 --bounded_du 0.15 \
    --w_sum 0.5 --w_shape 5.0 --w_U 0.5 --w_radius 5.0 \
    --progress_every 50
  echo EXIT=\$?
"

# balanced
OUT_SUB="$RAWDIR/v18_weight_balanced"
mkdir -p "$OUT_SUB"
launch_run "v18_weight_balanced" bash -c "
  cd '$OUT_SUB'
  '$PROJDIR/inc_ricci_v18' --nx 16 --ny 16 \
    --radius_seed 101 --target_dist lognormal --target_delta 0.25 \
    --target_radius_ratio 5 --coord_mode soft_relax --target_phi 0.84 \
    --relax_steps 20000 --use_voropp 1 --compute_nm_theta 1 \
    --decompose_theta 1 --apply_to_polydisperse 1 \
    --cycles 200 --shape_steps_per_cycle 50 --ricci_steps_per_cycle 50 \
    --shape_gamma 0.05 --bounded_du 0.15 \
    --w_sum 2.0 --w_shape 2.0 --w_U 1.0 --w_radius 5.0 \
    --progress_every 50
  echo EXIT=\$?
"

# weak_overlap
OUT_SUB="$RAWDIR/v18_weight_weak_overlap"
mkdir -p "$OUT_SUB"
launch_run "v18_weight_weak_overlap" bash -c "
  cd '$OUT_SUB'
  '$PROJDIR/inc_ricci_v18' --nx 16 --ny 16 \
    --radius_seed 101 --target_dist lognormal --target_delta 0.25 \
    --target_radius_ratio 5 --coord_mode soft_relax --target_phi 0.84 \
    --relax_steps 20000 --use_voropp 1 --compute_nm_theta 1 \
    --decompose_theta 1 --apply_to_polydisperse 1 \
    --cycles 200 --shape_steps_per_cycle 50 --ricci_steps_per_cycle 50 \
    --shape_gamma 0.05 --bounded_du 0.15 \
    --w_sum 1.0 --w_shape 2.0 --w_U 0.05 --w_radius 5.0 \
    --progress_every 50
  echo EXIT=\$?
"

# Wait for v18 scans to finish before starting v17 (to not overload)
echo "Waiting for v18 scans to complete..."
wait_all

# ===========================================================
# SECTION G: v17 Operator Atlas Extension
# ===========================================================
echo "=== SECTION G: v17 Operator Atlas ==="

# 16x16 basic operators (N=256)
OUT_SUB="$RAWDIR/v17_op_16_basic_more"
mkdir -p "$OUT_SUB"
launch_run "v17_operator_16_basic" bash -c "
  cd '$OUT_SUB'
  '$PROJDIR/inc_ricci_v17' --nx 16 --ny 16 \
    --seed 5000 --radius_seed 101 --same_radii_all_starts 1 \
    --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 \
    --coord_mode soft_relax --target_phi 0.84 --relax_steps 15000 \
    --use_voropp 1 --compute_nm_theta 1 --decompose_theta 1 \
    --operator_multistart 1 --operator_set basic \
    --threads 20 --operator_multistart_count 128 \
    --operator_multistart_out operator_multistart_v17_16_basic_more.csv \
    --progress_every 0
  echo EXIT=\$?
"

# 24x24 basic operators (N=576)
OUT_SUB="$RAWDIR/v17_op_24_basic_more"
mkdir -p "$OUT_SUB"
launch_run "v17_operator_24_basic" bash -c "
  cd '$OUT_SUB'
  '$PROJDIR/inc_ricci_v17' --nx 24 --ny 24 \
    --seed 6000 --radius_seed 202 --same_radii_all_starts 1 \
    --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 \
    --coord_mode soft_relax --target_phi 0.84 --relax_steps 25000 \
    --use_voropp 1 --compute_nm_theta 1 --decompose_theta 1 \
    --operator_multistart 1 --operator_set basic \
    --threads 20 --operator_multistart_count 48 \
    --operator_multistart_out operator_multistart_v17_24_basic_more.csv \
    --progress_every 0
  echo EXIT=\$?
"

# ===========================================================
# SECTION H: v17 Fine Compression Scan
# ===========================================================
echo "=== SECTION H: Fine Compression ==="
OUT_SUB="$RAWDIR/v17_fine_compression"
mkdir -p "$OUT_SUB"
launch_run "v17_fine_compression" bash -c "
  cd '$OUT_SUB'
  '$PROJDIR/inc_ricci_v17' --nx 16 --ny 16 \
    --radius_seed 101 --same_radii_all_starts 1 \
    --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 \
    --coord_mode soft_relax --relax_steps 20000 \
    --use_voropp 1 --compute_nm_theta 1 --decompose_theta 1 \
    --fine_compression_scan 1 \
    --fine_phi_list 0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92 \
    --fine_compression_seeds 1,2,3,4,5,6,7,8 \
    --fine_compression_out fine_compression_v17_16_fine.csv \
    --progress_every 0
  echo EXIT=\$?
"

echo "Waiting for v17 scans to complete..."
wait_all

echo "=== ALL SCANS COMPLETE ==="
echo "End time: $(date)"
