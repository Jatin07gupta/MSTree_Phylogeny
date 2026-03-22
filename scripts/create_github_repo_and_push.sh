#!/usr/bin/env bash
# Create https://github.com/Jatin07gupta/MSTree_Phylogeny (if missing) and push main.
#
# Why this exists: `git push` alone cannot create a GitHub repository or authenticate
# without credentials. This script uses the GitHub REST API + a Personal Access Token.
#
# Prerequisites:
#   1. GitHub account: Jatin07gupta
#   2. A Personal Access Token (classic) with scope "repo", OR a token that can
#      create repositories and push (see README). Fine-grained tokens are often
#      per-repository — create an empty repo on the website first if API create fails.
#
# Usage (from anywhere):
#   export GITHUB_TOKEN=ghp_your_token_here
#   bash /path/to/clnj_cpp/scripts/create_github_repo_and_push.sh
#
# Security: Never commit the token. Revoke the token after use if it was temporary.

set -euo pipefail

OWNER="Jatin07gupta"
REPO="MSTree_Phylogeny"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ -z "${GITHUB_TOKEN:-}" ]]; then
  echo "Error: GITHUB_TOKEN is not set."
  echo "Create a token: GitHub → Settings → Developer settings → Personal access tokens."
  echo "Classic token: enable scope \"repo\" (for private) or at least \"public_repo\" for public repos."
  exit 1
fi

API="https://api.github.com"
AUTH=(-H "Authorization: Bearer ${GITHUB_TOKEN}" -H "Accept: application/vnd.github+json" -H "X-GitHub-Api-Version: 2022-11-28")

echo "==> Checking local git repo: ${ROOT}"
if [[ ! -d "${ROOT}/.git" ]]; then
  echo "Error: No .git directory under ${ROOT}"
  exit 1
fi

code="$(curl -sS -o /tmp/gh_repo_check.json -w "%{http_code}" "${AUTH[@]}" "${API}/repos/${OWNER}/${REPO}")"

if [[ "${code}" == "200" ]]; then
  echo "==> Remote repository already exists: ${OWNER}/${REPO}"
elif [[ "${code}" == "404" ]]; then
  echo "==> Creating public repository ${OWNER}/${REPO} via API..."
  create_body='{"name":"'"${REPO}"'","description":"Multi-model distance-based phylogenetic pipeline (C++)","private":false,"has_issues":true,"has_wiki":false,"auto_init":false}'
  create_code="$(curl -sS -o /tmp/gh_repo_create.json -w "%{http_code}" -X POST "${AUTH[@]}" "${API}/user/repos" -d "${create_body}")"
  if [[ "${create_code}" != "201" ]]; then
    echo "Error: GitHub API returned HTTP ${create_code} when creating repo."
    cat /tmp/gh_repo_create.json 2>/dev/null || true
    echo ""
    echo "Common fixes:"
    echo "  - Use a classic PAT with \"repo\" scope."
    echo "  - Or create the empty repo manually on github.com, then re-run this script (it will only push)."
    exit 1
  fi
  echo "==> Repository created."
elif [[ "${code}" == "401" ]] || [[ "${code}" == "403" ]]; then
  echo "Error: GitHub API returned HTTP ${code} (bad token or missing scope)."
  cat /tmp/gh_repo_check.json 2>/dev/null || true
  exit 1
else
  echo "Error: Unexpected HTTP ${code} when checking repo."
  cat /tmp/gh_repo_check.json 2>/dev/null || true
  exit 1
fi

echo "==> Pushing branch main to origin (authenticated URL, not saved in config)..."
cd "${ROOT}"
GIT_TERMINAL_PROMPT=0 git remote remove origin 2>/dev/null || true
GIT_TERMINAL_PROMPT=0 git remote add origin "https://github.com/${OWNER}/${REPO}.git"

# Push using token in URL. Do NOT use --set-upstream with this URL — Git would store the
# token inside .git/config. Push anonymously to the auth URL, then track origin/main.
if ! GIT_TERMINAL_PROMPT=0 git push "https://oauth2:${GITHUB_TOKEN}@github.com/${OWNER}/${REPO}.git" main; then
  echo "Error: git push failed. Check token has \"repo\" (or contents: write) and repo name matches."
  exit 1
fi

GIT_TERMINAL_PROMPT=0 git fetch origin
GIT_TERMINAL_PROMPT=0 git branch --unset-upstream 2>/dev/null || true
GIT_TERMINAL_PROMPT=0 git branch -u "origin/main" main

echo "==> Done. Open: https://github.com/${OWNER}/${REPO}"
echo "==> Remote 'origin' is set to https (without token). Configure credential helper for future pushes."
