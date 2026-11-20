#!/usr/bin/env bash
# Publish the MyST report (index.md) to GitHub Pages as static HTML.
#
# Builds locally (the site needs the annexed figures, which a CI runner would
# have to fetch from LFS) and force-pushes the result to the `gh-pages` branch,
# a deploy-only branch that holds no history worth keeping.
#
# The HDF5 output archives are stripped from the site: one exceeds GitHub's
# 100 MB file limit, and all of them are available from the repository via
# git-annex (see README "Getting the data"). Their links are pointed there.
#
# Usage: tools/publish_pages.sh            (from the project root)
set -euo pipefail

REPO_URL="git@github.com:dwoodward-0/LZLightCone.git"
BASE_URL="/LZLightCone"
DATA_URL="https://github.com/dwoodward-0/LZLightCone#getting-the-data"

[ -f astra.yaml ] || { echo "run from the project root" >&2; exit 1; }

rm -rf _build/html
BASE_URL="$BASE_URL" myst build --html

site=_build/html
# Strip HDF5 archives and redirect their links.
find "$site" -type f -name '*.hdf5' -delete
grep -rlE "${BASE_URL}/build/[A-Za-z0-9_]+-[0-9a-f]{32}\.hdf5" "$site" \
  | while read -r f; do
      sed -i '' -E "s|${BASE_URL}/build/[A-Za-z0-9_]+-[0-9a-f]{32}\.hdf5|${DATA_URL}|g" "$f"
    done
touch "$site/.nojekyll"
# Hydration fix: myst 1.13 bakes the page's route as "/index" without the base
# path. The Remix client compares it with window.location.pathname
# ("/LZLightCone/") and reloads on any mismatch, which loops forever. Bake in
# the served path instead (single-page site: index.html only).
sed -i '' "s|__remixContext = {\"url\":\"/index\"|__remixContext = {\"url\":\"${BASE_URL}/\"|" "$site/index.html"
grep -q "__remixContext = {\"url\":\"${BASE_URL}/\"" "$site/index.html" \
  || { echo "refusing to publish: hydration URL patch did not apply" >&2; exit 1; }
# Stable paths for the README's embedded headline figures.
mkdir -p "$site/figures"
cp results/baseline/sensitivity_vs_lz_validation.png \
   results/baseline/plr_fc_exposure_comparison.png "$site/figures/"

# Guards: nothing withheld from the public repo, nothing over GitHub's limit.
if find "$site" -type f \( -iname '*.pdf' -o -iname '*.py' -o -iname '*.hdf5' \) | grep -q .; then
  echo "refusing to publish: PDF/.py/HDF5 found in site" >&2; exit 1
fi
if find "$site" -type f -size +95M | grep -q .; then
  echo "refusing to publish: file over 95 MB in site" >&2; exit 1
fi

src_sha=$(git rev-parse --short HEAD)
# Push from a throwaway repo so the project's git-annex filters never apply.
tmp=$(mktemp -d)
cp -R "$site"/. "$tmp"/
(
  cd "$tmp"
  git init -q -b gh-pages
  git add -A
  git commit -q -m "Publish report (source ${src_sha})"
  git push -q -f "$REPO_URL" gh-pages
)
rm -rf "$tmp"
echo "published: https://dwoodward-0.github.io/LZLightCone/"
