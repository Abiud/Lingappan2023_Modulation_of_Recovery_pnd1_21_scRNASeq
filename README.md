
# Modulation of Recovery from Neonatal Hyperoxic Lung Injury by Sex as a Biological Variable
Please refer to the GEO dataset for a link to the publication.

The book folder contains quarto files with a step by step analysis.
The R folder conatins all of the functions used on each step. 

## GEO
[GSE237944](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237944)

## Docker
A docker file is included. This image installs all the necessary packages and code to reproduce our results. 

## Run Locally
R 4.2.3 and the package renv is recommended so the user can replicate the resulsts in the publication.

Clone the project, download processed files from the GEO into their respective folders. Each raw_data subfolder should have these 4 files:
* barcodes.tsv.gz
* features.tsv.gz
* matrix.mtx.gz
* metadata.csv
Rename these files to match the above names.

Open an R session and make the project directory your working directory.
```R
  setwd("path/to/this/Lingappan2023_Modulation_of_Recovery_pnd1_21")
```

Install needed packages to run pipeline.
```R
  renv::restore()
```

## Shiny Apps

[PND 1 & PND 21 Integrated Seurat objects - ShinyCell](https://abiudcantu.shinyapps.io/pnd1_21_seurat/)

More Shiny Apps to visualize our data are in our [website!](https://www.lingappanlab.com/)

## Questions & Comments

Please contact Abiud (the owner of this repo) for any questions or comments.
