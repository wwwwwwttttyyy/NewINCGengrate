#!/bin/bash
# Quick check on gamma scan results
for g in 0.02 0.05 0.10 0.25 0.50 1.0; do
  dir="/mnt/d/Codes/airesearch/NewINCGengrate/old_program_exploration/raw/v18_poly_gamma_${g}"
  if [ -f "$dir/coupled_cycle_trajectory.csv" ]; then
    lines=$(wc -l < "$dir/coupled_cycle_trajectory.csv")
    # Count accepted cycles (column 15)
    accepted=$(awk -F, 'NR>1 && $15==1' "$dir/coupled_cycle_trajectory.csv" | wc -l)
    # Get first and last values
    first_sum=$(awk -F, 'NR==2{print $4}' "$dir/coupled_cycle_trajectory.csv")
    last_sum=$(awk -F, 'END{print $4}' "$dir/coupled_cycle_trajectory.csv")
    first_sh=$(awk -F, 'NR==2{print $5}' "$dir/coupled_cycle_trajectory.csv")
    last_sh=$(awk -F, 'END{print $5}' "$dir/coupled_cycle_trajectory.csv")
    echo "gamma=$g cycles=$((lines-1)) accepted=$accepted sum: $first_sum -> $last_sum shape: $first_sh -> $last_sh"
  else
    echo "gamma=$g NOT_READY"
  fi
done
echo "=== Bound scan ==="
for b in 0.03 0.05 0.15 0.25 0.40; do
  dir="/mnt/d/Codes/airesearch/NewINCGengrate/old_program_exploration/raw/v18_poly_bound_${b}"
  if [ -f "$dir/coupled_cycle_trajectory.csv" ]; then
    lines=$(wc -l < "$dir/coupled_cycle_trajectory.csv")
    accepted=$(awk -F, 'NR>1 && $15==1' "$dir/coupled_cycle_trajectory.csv" | wc -l)
    echo "bound_du=$b cycles=$((lines-1)) accepted=$accepted"
  else
    echo "bound_du=$b NOT_READY"
  fi
done
echo "=== Weight scan ==="
for w in shape_biased balanced weak_overlap; do
  dir="/mnt/d/Codes/airesearch/NewINCGengrate/old_program_exploration/raw/v18_weight_${w}"
  if [ -f "$dir/coupled_cycle_trajectory.csv" ]; then
    lines=$(wc -l < "$dir/coupled_cycle_trajectory.csv")
    accepted=$(awk -F, 'NR>1 && $15==1' "$dir/coupled_cycle_trajectory.csv" | wc -l)
    echo "weight=$w cycles=$((lines-1)) accepted=$accepted"
  else
    echo "weight=$w NOT_READY"
  fi
done
