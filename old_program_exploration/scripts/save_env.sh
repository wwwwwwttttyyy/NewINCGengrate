#!/bin/bash
OUT=/mnt/d/Codes/airesearch/NewINCGengrate/old_program_exploration/reports/environment.txt
{
echo "=== Environment Info ==="
echo "date: $(date)"
echo "pwd: $(pwd)"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate analyse
echo "conda env: $CONDA_DEFAULT_ENV"
g++ --version | head -1
echo "voro++ path: $(which voro++)"
echo
echo "=== File sizes ==="
ls -lh /mnt/d/Codes/airesearch/NewINCGengrate/inc_ricci_v16
ls -lh /mnt/d/Codes/airesearch/NewINCGengrate/inc_ricci_v17
ls -lh /mnt/d/Codes/airesearch/NewINCGengrate/inc_ricci_v18
ls -lh /mnt/d/Codes/airesearch/NewINCGengrate/inc_ricci_v16.cpp
ls -lh /mnt/d/Codes/airesearch/NewINCGengrate/inc_ricci_v17.cpp
ls -lh /mnt/d/Codes/airesearch/NewINCGengrate/inc_ricci_v18.cpp
} > "$OUT" 2>&1
cat "$OUT"