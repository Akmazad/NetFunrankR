# NetFunrankR dashboard container.
# Build from the repository root:
#   docker build -f deploy/Dockerfile -t netfunrankr .
#   docker run -p 7860:7860 netfunrankr          # then open http://localhost:7860
# Works unchanged on Hugging Face Spaces (Docker SDK, port 7860), Google Cloud
# Run, Fly.io, Render or any Docker host. Configure with the NFR_* variables
# documented in deploy/app/app.R.
FROM rocker/r-ver:4.4.2

RUN apt-get update && apt-get install -y --no-install-recommends \
        libglpk-dev libxml2-dev libcurl4-openssl-dev libssl-dev \
        libpng-dev zlib1g-dev pandoc \
    && rm -rf /var/lib/apt/lists/*

# Dependencies (binary packages from the Posit Package Manager configured in
# rocker images, Bioconductor packages through BiocManager)
RUN Rscript -e 'install.packages(c("BiocManager", "igraph", "Matrix", \
        "Rcpp", "rlang", "ggplot2", "shiny", "bslib", "bsicons", "plotly", \
        "reactable", "visNetwork", "shinyvalidate", "htmltools", "mirai", \
        "promises", "rmarkdown", "knitr"))' \
 && Rscript -e 'BiocManager::install(c("S4Vectors", "BiocParallel"), \
        update = FALSE, ask = FALSE)'

# The package itself (build context: repository root; see .dockerignore)
COPY DESCRIPTION NAMESPACE LICENSE /src/NetFunrankR/
COPY R /src/NetFunrankR/R
COPY src /src/NetFunrankR/src
COPY data /src/NetFunrankR/data
COPY inst /src/NetFunrankR/inst
COPY man /src/NetFunrankR/man
RUN R CMD INSTALL --no-test-load /src/NetFunrankR && rm -rf /src

COPY deploy/app /app

# Hugging Face Spaces run containers as user 1000
RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user NFR_WORKERS=2 PORT=7860
EXPOSE 7860
CMD ["Rscript", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = as.integer(Sys.getenv('PORT', '7860')), launch.browser = FALSE)"]
