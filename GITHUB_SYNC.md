# GitHub Sync Notes

This repository is intended to sync the reproducible research code, reports, and scripts between machines.

## What Goes In Git

Tracked:

- `inc_ricci*.cpp`
- `REPORT*.md`
- analysis scripts and reports under `old_program_exploration/`
- small hand-written support files

Ignored:

- compiled binaries such as `inc_ricci_v18`, `*.exe`, and ASAN builds
- generated run directories such as `out_*`, `coord_*`, `theta_*`, `voro_*`
- root-level generated CSV/TXT trajectory and summary files

Large generated outputs should be regenerated from commands in the reports, archived separately, or moved to Git LFS later if needed.

## First Push From This Machine

Create an empty GitHub repository in the browser, then run:

```bash
git remote add origin https://github.com/<USER>/<REPO>.git
git branch -M main
git push -u origin main
```

If GitHub asks for credentials over HTTPS, use a GitHub personal access token instead of your password.

## On The Other Machine

```bash
git clone https://github.com/<USER>/<REPO>.git
cd <REPO>
```

Build under WSL conda `analyse`:

```bash
g++ -O3 -std=c++17 inc_ricci_v18.cpp -o inc_ricci_v18 \
  -I/home/tianyu/miniconda3/envs/analyse/include \
  -I/home/tianyu/miniconda3/envs/analyse/include/voro++ \
  -L/home/tianyu/miniconda3/envs/analyse/lib \
  -lvoro++ \
  -Wl,-rpath,/home/tianyu/miniconda3/envs/analyse/lib \
  -pthread
```

## Daily Two-Machine Workflow

Before starting work:

```bash
git pull
```

After editing code or reports:

```bash
git status
git add inc_ricci_v*.cpp REPORT_V*.md old_program_exploration GITHUB_SYNC.md .gitignore
git commit -m "Describe the change"
git push
```

Do not commit generated output directories unless there is a specific reason.

## One-Command Sync Scripts

PowerShell on Windows:

```powershell
.\scripts\git_update.ps1
# run experiments or edit code
.\scripts\git_upload.ps1 -Message "describe the update"
```

WSL / bash:

```bash
bash scripts/git_update.sh
# run experiments or edit code
bash scripts/git_upload.sh "describe the update"
```

`git_update` fetches and rebases onto `origin/main` with `--autostash`.
`git_upload` stages tracked and unignored changes, creates a commit if needed,
rebases onto the latest `origin/main`, and pushes.

The scripts do not force-push. If both machines changed the same tracked file,
Git may stop for a conflict; resolve it manually, then continue with
`git rebase --continue` and rerun the upload script.
