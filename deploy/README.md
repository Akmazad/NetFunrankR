# Deploying the NetFunrankR dashboard

The dashboard is the function `netFunrankRApp()` of the package, so any host
that can run an R Shiny app can serve it. Everything a deployment needs is in
this folder:

| File | Purpose |
|---|---|
| `app/app.R` | entry point; configured with `NFR_*` environment variables |
| `Dockerfile` | container image (Hugging Face Spaces, Cloud Run, any Docker host) |
| `make_hf_space.sh` | assembles a ready-to-push Hugging Face Space |
| `deploy_posit.R` | deploys to shinyapps.io / writes the manifest for Posit Connect Cloud |

## Resource needs

Measured on the bundled STRING network (16,024 proteins):

| Libraries offered | Peak memory per R process |
|---|---|
| KEGG or WikiPathways | about 0.5 GB |
| Reactome (1,700 sets) | about 0.65 GB |

With background workers (`NFR_WORKERS >= 1`, the default), the app process
and each worker hold their own copy of the network, so plan for
(1 + workers) x the figure above. On servers with 1 GB of memory, set
`NFR_ASYNC=false` (compute in the app process) and/or offer only KEGG and
WikiPathways (`NFR_LIBRARIES="KEGG 2021,WikiPathways 2024"`).

## Free hosting options

Free tiers change; check the provider's current limits before relying on
them.

### 1. Hugging Face Spaces (recommended)

Docker Spaces on the free CPU hardware have ample memory for all libraries
and keep a public URL. Spaces go to sleep after a period without visitors
and wake on the next visit (first load then takes a minute).

```bash
# from the repository root
bash deploy/make_hf_space.sh <hf-user>/netfunrankr hf-space
# create an empty Space (SDK: Docker) at https://huggingface.co/new-space
cd hf-space
git init && git remote add origin https://huggingface.co/spaces/<hf-user>/netfunrankr
git add . && git commit -m "NetFunrankR dashboard" && git push -u origin main
```

The Space builds the image from `Dockerfile` (about 10-15 minutes the first
time) and serves the app at `https://huggingface.co/spaces/<hf-user>/netfunrankr`.

### 2. Posit Connect Cloud

Deploys directly from a public GitHub repository.

```r
remotes::install_github("Akmazad/NetFunrankR")   # so the source is recorded
source("deploy/deploy_posit.R")                   # writes deploy/app/manifest.json
```

Commit `deploy/app/manifest.json`, sign in at https://connect.posit.cloud
with GitHub, choose *Publish > Shiny*, select the repository and
`deploy/app/app.R`.

### 3. shinyapps.io

The free plan limits active hours per month and uses small instances, so
set `NFR_ASYNC=false` there.

```r
remotes::install_github("Akmazad/NetFunrankR")
rsconnect::setAccountInfo(name = "<account>", token = "<token>", secret = "<secret>")
```

```bash
Rscript deploy/deploy_posit.R shinyapps <account>
```

### 4. Your own server or institutional VM

```bash
docker build -f deploy/Dockerfile -t netfunrankr .
docker run -d -p 80:7860 -e NFR_WORKERS=4 netfunrankr
```

or install the package on a machine running the free Shiny Server and copy
`deploy/app/` into its site directory.

## Customising a deployment

Without changing code:

```bash
docker run -p 7860:7860 \
  -e NFR_LIBRARIES="KEGG 2021" \
  -e NFR_GMT_DIR=/data/gmt -e NFR_EDGE_DIR=/data/networks \
  -v /path/to/my/resources:/data netfunrankr
```

Every `.gmt` file in `NFR_GMT_DIR` becomes an extra library and every
two-column `.csv`/`.tsv` edge list in `NFR_EDGE_DIR` an extra interactome
(any organism or identifier type). For full control, write your own
`app.R` calling `netFunrankRApp(geneSets = ..., networks = ...,
exampleGenes = ..., geneLabel = ...)`; see `?runNetFunrankRApp`.

## Not currently possible

Fully static, server-less hosting (Shinylive/webR on GitHub Pages) is not
supported, because the package has compiled C++ code and Bioconductor
dependencies that are not available as WebAssembly builds.
