#!/usr/bin/env bash
set -euo pipefail

MESSAGE="${1:-}"
REMOTE="${REMOTE:-origin}"
BRANCH="${BRANCH:-main}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

git rev-parse --is-inside-work-tree >/dev/null

if [[ -z "${MESSAGE}" ]]; then
  MESSAGE="sync update $(date '+%Y-%m-%d %H:%M:%S')"
fi

echo "Repository: ${REPO_ROOT}"
echo "Staging tracked and unignored changes ..."
git add -A

if ! git diff --cached --quiet; then
  echo "Creating commit: ${MESSAGE}"
  git commit -m "${MESSAGE}"
else
  echo "No staged changes to commit."
fi

echo "Rebasing on latest ${REMOTE}/${BRANCH} before push ..."
git fetch "${REMOTE}"
git pull --rebase --autostash "${REMOTE}" "${BRANCH}"

echo "Pushing to ${REMOTE}/${BRANCH} ..."
git push "${REMOTE}" "${BRANCH}"

echo
echo "Upload complete."
git status -sb
