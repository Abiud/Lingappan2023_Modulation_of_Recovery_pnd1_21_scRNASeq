FROM rocker/r-ver:4.2.3@sha256:7c89f3b30b3f998215a74bfd992af90346f34b4f3d233cdfb0c93634f3b9ecf4

LABEL org.opencontainers.image.licenses="GPL-2.0-or-later" \
    org.opencontainers.image.source="https://github.com/rocker-org/rocker-versioned2" \
    org.opencontainers.image.vendor="Rocker Project" \
    org.opencontainers.image.authors="Carl Boettiger <cboettig@ropensci.org>"

ENV S6_VERSION=v2.1.0.2
ENV RSTUDIO_VERSION=2023.03.0+386
ENV DEFAULT_USER=rstudio
ENV PANDOC_VERSION=default
ENV QUARTO_VERSION=default
ENV RENV_VERSION 1.0.0
ENV RENV_PATHS_LIBRARY renv/library

RUN /rocker_scripts/install_rstudio.sh
RUN /rocker_scripts/install_pandoc.sh
RUN /rocker_scripts/install_quarto.sh
RUN apt-get update && apt-get install -y \
    libudunits2-dev \
    libcairo2-dev \
    libxt-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    cmake \
    libharfbuzz-dev \
    libfribidi-dev \
    && rm -rf /var/lib/apt/lists/*
RUN R -e "install.packages('renv')"

ADD renv.lock /home/rstudio/
ADD book /home/rstudio/book
ADD data /home/rstudio/data
ADD R /home/rstudio/R
ADD raw_data /home/rstudio/raw_data
ADD renv /home/rstudio/renv

# RUN R -e "devtools::install_version('speedglm', '0.3-4', repos = 'https://packagemanager.rstudio.com/cran/2023-03-31')"
RUN R -e "renv::restore()"
# RUN R -e "remotes::install_github('satijalab/seurat-wrappers')"

EXPOSE 8787

CMD ["/init"]
