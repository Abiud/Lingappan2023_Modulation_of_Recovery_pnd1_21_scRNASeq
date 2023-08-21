# Modulation of Recovery from Neonatal Hyperoxic Lung Injury by Sex as a Biological Variable

Please refer to the GEO dataset for a link to the publication.

The book folder contains quarto files with a step by step analysis.
The R folder conatins all of the functions used on each step.

## GEO

[GSE237944](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237944)

## Docker

A docker file is included. This image installs all the necessary packages and code to reproduce our results. The image is also uploaded to [Docker hub](https://hub.docker.com/repository/docker/abiud/lingappan2023_modulation_of_recovery/general)

## Run Locally

R 4.2.3 and the package renv is recommended so the user can replicate the resulsts in the publication.

Clone the project, download processed files from the GEO into their respective folders. Each raw_data subfolder should have these 4 files:

- barcodes.tsv.gz
- features.tsv.gz
- matrix.mtx.gz
- metadata.csv
  Rename these files to match the above names.

Install needed packages to run pipeline.

```R
  renv::restore()
```

## Shiny Apps

[PND 1 & PND 21 Integrated Seurat objects - ShinyCell](https://abiudcantu.shinyapps.io/pnd1_21_seurat/)

[PND 1 CellChat](https://abiudcantu.shinyapps.io/pnd1_immune_endo/)

[PND 21 CellChat](https://abiudcantu.shinyapps.io/pnd21_immune_endo_cellchat/)

More Shiny Apps to visualize our data are in our [website!](https://www.lingappanlab.com/)

## Questions & Comments

Please contact Abiud (the owner of this repo) for any questions or comments.
