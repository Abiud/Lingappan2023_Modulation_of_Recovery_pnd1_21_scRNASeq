#' Read 10X data into a Seurat object
#'
#' @param sampleDir A character string specifying the directory containing the sample data
#' @param projectName A character string specifying the name of the project
#'
#' @return A Seurat object containing the gene expression data
#'
#' @export
Read10XToSeurat <- function(sample) {
    tmpData <- Seurat::Read10X(data.dir = paste0("./raw_data/", sample, "/"))
    seuratSample <- Seurat::CreateSeuratObject(counts = tmpData, project = sample)
    seuratSample$percent.mt <- Seurat::PercentageFeatureSet(seuratSample, pattern = "^mt-")
    seuratSample <- Seurat::SCTransform(seuratSample, vars.to.regress = "percent.mt", verbose = FALSE)
    dir.create(file.path("./figures", age, "QC", sample), recursive = TRUE, showWarnings = FALSE)
    seuratSample <- GetSeuratClusters(seuratSample, stepName = "10XToSeurat")
    return(seuratSample)
}
