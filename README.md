# INC-Ricci Research Workspace

Single-file C++ research programs and reports for the INC/Ricci topology experiments.

Primary current files:

- `inc_ricci_v18.cpp`: validated operator-block version
- `inc_ricci_v19.cpp`: work-in-progress failure-forensics version
- `REPORT_V18.md`: latest completed implementation report
- `old_program_exploration/reports/old_program_exploration_report.md`: recent failure analysis of v17/v18 runs
- `GITHUB_SYNC.md`: two-machine Git/GitHub workflow

The repository intentionally ignores compiled binaries and generated run outputs. Rebuild executables from source and regenerate large outputs from the report commands.

WSL/voro++ build example:

```bash
g++ -O3 -std=c++17 inc_ricci_v18.cpp -o inc_ricci_v18 \
  -I/home/tianyu/miniconda3/envs/analyse/include \
  -I/home/tianyu/miniconda3/envs/analyse/include/voro++ \
  -L/home/tianyu/miniconda3/envs/analyse/lib \
  -lvoro++ \
  -Wl,-rpath,/home/tianyu/miniconda3/envs/analyse/lib \
  -pthread
```
