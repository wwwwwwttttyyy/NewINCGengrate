#!/bin/bash
# Run v18 gamma scan - properly handles CWD output
PROJDIR="/mnt/d/Codes/airesearch/NewINCGengrate"
RAWDIR="$PROJDIR/old_program_exploration/raw"
LOGDIR="$PROJDIR/old_program_exploration/logs"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate analyse

cd "$PROJDIR"

for gamma in 0 0.02 0.05 0.10 0.25 0.50 1.0; do
  echo "--- gamma=$gamma START $(date) ---"
  OUTDIR="$RAWDIR/v18_poly_gamma_${gamma}"
  mkdir -p "$OUTDIR"
  cd "$OUTDIR"
  "$PROJDIR/inc_ricci_v18" --nx 16 --ny 16 --radius_seed 101 --target_dist lognormal --target_delta 0.25 --target_radius_ratio 5 --coord_mode soft_relax --target_phi 0.84 --relax_steps 20000 --use_voropp 1 --compute_nm_theta 1 --decompose_theta 1 --apply_to_polydisperse 1 --cycles 200 --shape_steps_per_cycle 50 --ricci_steps_per_cycle 50 --shape_gamma ${gamma} --bounded_du 0.10 --progress_every 50 > "$LOGDIR/v18_poly_gamma_${gamma}.log" 2>&1
  echo "--- gamma=$gamma DONE exit=$? $(date) ---"
  cd "$PROJDIR"
done

echo "=== GAMMA SCAN COMPLETE ==="