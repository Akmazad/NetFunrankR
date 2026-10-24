#!/usr/bin/env bash
# Assemble a ready-to-push Hugging Face Space (Docker SDK) for the dashboard.
#
# Usage (from the repository root):
#   bash deploy/make_hf_space.sh <hf-username>/<space-name> [output-dir]
# then
#   cd <output-dir> && git init && git remote add origin \
#     https://huggingface.co/spaces/<hf-username>/<space-name> && \
#     git add . && git commit -m "NetFunrankR dashboard" && git push -u origin main
set -euo pipefail
SPACE="${1:?usage: make_hf_space.sh <user>/<space-name> [output-dir]}"
OUT="${2:-hf-space}"
rm -rf "$OUT" && mkdir -p "$OUT/deploy"
cp -r DESCRIPTION NAMESPACE LICENSE R src data inst man "$OUT"/
rm -f "$OUT"/src/*.o "$OUT"/src/*.so "$OUT"/src/*.dll
cp -r deploy/app "$OUT/deploy/app"
cp deploy/Dockerfile "$OUT/Dockerfile"
cp .dockerignore "$OUT/.dockerignore"
cat > "$OUT/README.md" <<EOF
---
title: NetFunrankR
emoji: 🧬
colorFrom: blue
colorTo: green
sdk: docker
app_port: 7860
pinned: false
license: mit
short_description: Decomposed network-proximity enrichment of gene lists
---

# NetFunrankR dashboard

Interactive version of the R/Bioconductor package
[NetFunrankR](https://github.com/Akmazad/NetFunrankR): test whether a gene
list is close to pathways in the human interactome, with the evidence split
into annotation overlap and network (linker) genes.

Space: https://huggingface.co/spaces/${SPACE}
EOF
echo "Space files written to $OUT"
