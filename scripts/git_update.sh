#!/usr/bin/env bash
set -euo pipefail

REMOTE="${1:-origin}"
BRANCH="${2:-main}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

git rev-parse --is-inside-work-tree >/dev/null

echo "Repository: ${REPO_ROOT}"
echo "Updating from ${REMOTE}/${BRANCH} ..."

git fetch "${REMOTE}"
git pull --rebase --autostash "${REMOTE}" "${BRANCH}"

echo
echo "Update complete."
git status -sb
