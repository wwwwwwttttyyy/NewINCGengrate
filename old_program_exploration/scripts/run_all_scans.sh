#!/bin/bash
# Old Program Exploration - Master Run Script
# v18 writes output to CWD, so we cd to each output dir before running
set -e
source ~/miniconda3/etc/profile.d/conda.sh
conda activate analyse
PROJDIR="/mnt/d/Codes/airesearch/NewINCGengrate"
OUTDIR="$PROJDIR/old_program_exploration"
LOGDIR="$OUTDIR/logs"
RAWDIR="$OUTDIR/raw"
cd "$PROJDIR"

echo "=== Environment ==="
echo "date: $(date)"
echo "pwd: $(pwd)"
echo "conda env: $CONDA_DEFAULT_ENV"
g++ --version | head -1
which voro++
echo "=== End Environment ==="

# SECTION C: v18 Gamma Scan
echo "=== SECTION C: v18 Gamma Scan ==="
for gamma in 0 0.02 0.05 0.10 0.25 0.50 1.0; do
  echo "--- gamma=$gamma START $(date) ---"
  OUT_SUB="$RAWDIR/v18_poly_gamma_${gamma}"
  mkdir -p "$OUT_SUB"
  cd "$OUT_SUB"
  $PROJDIR/inc_ricci_v18 --nx 16 --ny 16 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --coord_mode soft_relax --target_phi 0.84 --relax_steps 20000 --use_voropp 1 --compute_nm_theta 1 --decompose_theta 1 --apply_to_polydisperse 1 --cycles 200 --shape_steps_per_cycle 50 --ricci_steps_per_cycle 50 --shape_gamma ${gamma} --bounded_du 0.10 --progress_every 50 > "$LOGDIR/v18_poly_gamma_${gamma}.log" 2>&1
  echo "--- gamma=$gamma DONE exit=$? $(date) ---"
  cd "$PROJDIR"
done

# SECTION D: v18 Bounded_du Scan
echo "=== SECTION D: v18 Bounded_du Scan ==="
for bdu in 0.03 0.05 0.10 0.15 0.25 0.40; do
  echo "--- bounded_du=$bdu START $(date) ---"
  OUT_SUB="$RAWDIR/v18_poly_bound_${bdu}"
  mkdir -p "$OUT_SUB"
  cd "$OUT_SUB"
  $PROJDIR/inc_ricci_v18 --nx 16 --ny 16 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --coord_mode soft_relax --target_phi 0.84 --relax_steps 20000 --use_voropp 1 --compute_nm_theta 1 --decompose_theta 1 --apply_to_polydisperse 1 --cycles 200 --shape_steps_per_cycle 20 --ricci_steps_per_cycle 100 --shape_gamma 0.05 --bounded_du ${bdu} --progress_every 50 > "$LOGDIR/v18_poly_bound_${bdu}.log" 2>&1
  echo "--- bdu=$bdu DONE exit=$? $(date) ---"
  cd "$PROJDIR"
done

# SECTION E: v18 Weight Scan
echo "=== SECTION E: v18 Weight Scan ==="
for wscheme in shape_biased balanced weak_overlap; do
  echo "--- weight=$wscheme START $(date) ---"
  OUT_SUB="$RAWDIR/v18_weight_${wscheme}"
  mkdir -p "$OUT_SUB"
  cd "$OUT_SUB"
  WFLAGS=""
  case $wscheme in
    shape_biased) WFLAGS="--w_sum 0.5 --w_shape 5.0 --w_u 0.5 --w_radius 5.0" ;;
    balanced) WFLAGS="--w_sum 2.0 --w_shape 2.0 --w_u 1.0 --w_radius 5.0" ;;
    weak_overlap) WFLAGS="--w_sum 1.0 --w_shape 2.0 --w_u 0.05 --w_radius 5.0" ;;
  esac
  $PROJDIR/inc_ricci_v18 --nx 16 --ny 16 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --coord_mode soft_relax --target_phi 0.84 --relax_steps 20000 --use_voropp 1 --compute_nm_theta 1 --decompose_theta 1 --apply_to_polydisperse 1 --cycles 200 --shape_steps_per_cycle 50 --ricci_steps_per_cycle 50 --shape_gamma 0.05 --bounded_du 0.15 $WFLAGS --progress_every 50 > "$LOGDIR/v18_weight_${wscheme}.log" 2>&1
  echo "--- weight=$wscheme DONE exit=$? $(date) ---"
  cd "$PROJDIR"
done

# SECTION G: v17 Operator Atlas Extension
echo "=== SECTION G: v17 Operator Atlas ==="
echo "--- 16x16 basic START $(date) ---"
mkdir -p "$RAWDIR/v17_op_16_basic_more"
cd "$RAWDIR/v17_op_16_basic_more"
$PROJDIR/inc_ricci_v17 --nx 16 --ny 16 --seed 5000 --radius_seed 101 --same_radii_all_starts 1 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --coord_mode soft_relax --target_phi 0.84 --relax_steps 15000 --use_voropp 1 --compute_nm_theta 1 --decompose_theta 1 --operator_multistart 1 --operator_set basic --threads 20 --operator_multistart_count 128 --operator_multistart_out operator_multistart_v17_16_basic_more.csv --progress_every 0 > "$LOGDIR/v17_operator_16_basic.log" 2>&1
echo "--- 16x16 DONE exit=$? $(date) ---"
cd "$PROJDIR"

echo "--- 24x24 basic START $(date) ---"
mkdir -p "$RAWDIR/v17_op_24_basic_more"
cd "$RAWDIR/v17_op_24_basic_more"
$PROJDIR/inc_ricci_v17 --nx 24 --ny 24 --seed 6000 --radius_seed 202 --same_radii_all_starts 1 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --coord_mode soft_relax --target_phi 0.84 --relax_steps 25000 --use_voropp 1 --compute_nm_theta 1 --decompose_theta 1 --operator_multistart 1 --operator_set basic --threads 20 --operator_multistart_count 48 --operator_multistart_out operator_multistart_v17_24_basic_more.csv --progress_every 0 > "$LOGDIR/v17_operator_24_basic.log" 2>&1
echo "--- 24x24 DONE exit=$? $(date) ---"
cd "$PROJDIR"

# SECTION H: v17 Fine Compression Scan
echo "=== SECTION H: Fine Compression ==="
echo "--- fine compression START $(date) ---"
mkdir -p "$RAWDIR/v17_fine_compression"
cd "$RAWDIR/v17_fine_compression"
$PROJDIR/inc_ricci_v17 --nx 16 --ny 16 --radius_seed 101 --same_radii_all_starts 1 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --coord_mode soft_relax --relax_steps 20000 --use_voropp 1 --compute_nm_theta 1 --decompose_theta 1 --fine_compression_scan 1 --fine_phi_list 0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92 --fine_compression_seeds 1,2,3,4,5,6,7,8 --fine_compression_out fine_compression_v17_16_fine.csv --progress_every 0 > "$LOGDIR/v17_fine_compression.log" 2>&1
echo "--- fine compression DONE exit=$? $(date) ---"
cd "$PROJDIR"

echo "=== ALL SCANS COMPLETE ==="
echo "End time: $(date)"