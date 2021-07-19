# Docker inheritance
FROM bioconductor/bioconductor_docker:devel

# Update apt-get
RUN apt-get update \
## Remove packages in '/var/cache/' and 'var/lib'
## to remove side-effects of apt-get update
&& apt-get clean \
&& rm -rf /var/lib/apt/lists/*
# Install required Bioconductor packages mia and miaviz
RUN R -e 'BiocManager::install(c("mia","miaViz"))'