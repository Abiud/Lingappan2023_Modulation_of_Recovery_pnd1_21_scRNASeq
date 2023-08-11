#' Remove contamination from a Seurat object
#'
#' This function removes cells from a Seurat object that do not pass a set of quality control metrics, as well as a set of specified genes.
#'
#' @param seuratObj A Seurat object.
#' @param min_nUMI The lower threshold for the number of UMIs per cell (default is 500).
#' @param min_nGene The lower threshold for the number of genes per cell (default is 250).
#' @param min_log10GenesPerUMI The lower threshold for the log10(genes per UMI) metric (default is 0.80).
#' @param max_mitoRatio The upper threshold for the mitochondrial gene ratio metric (default is 0.05).
#' @param genesToRemove A character vector of gene names to remove (default is an empty vector).
#'
#' @return A preprocessed Seurat object.
#'
#' @export
RemoveContamination <- function(seuratObj, min_nUMI = 500, min_nGene = 250, min_log10GenesPerUMI = 0.80, max_mitoRatio = 0.05, genesToRemove = c()) {
    seuratObj <- RemoveGenes(seuratObj, genesToRemove)

    seuratObj <- FilterContamination(seuratObj, nUMI_l = min_nUMI, nGene_l = min_nGene, log10GenesPerUMI_l = min_log10GenesPerUMI, mitoRatio_u = max_mitoRatio)
    seuratObj <- Seurat::SCTransform(seuratObj, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = FALSE)
    seuratObj <- GetSeuratClusters(seuratObj, stepName = "removeCont")
    chosenResolution <- 0.6
    Idents(object = seuratObj) <- paste0("SCT_snn_res.", chosenResolution) # "SCT_snn_res.0.6"

    return(seuratObj)
}

#' Remove contaminated cells from a Seurat object
#'
#' @param seuratObj A Seurat object containing gene expression data for a specific sample
#' @param nUMI_l An integer specifying the lower threshold for UMI counts
#' @param nGene_l An integer specifying the lower threshold for gene counts
#' @param log10GenesPerUMI_l A numeric value specifying the lower threshold for log10(genes per UMI)
#' @param mitoRatio_u A numeric value specifying the upper threshold for mitochondrial gene expression
#' @param verbose A logical value specifying whether to print progress messages
#'
#' @return A new Seurat object with contaminated cells removed
#'
#' @export
FilterContamination <- function(seuratObj, nUMI_l = 500, nGene_l = 250, log10GenesPerUMI_l = 0.80, mitoRatio_u = 0.05) {
    sMetadata <- seuratObj@meta.data
    sMetadata$cells <- rownames(sMetadata)

    outlierTreshold_nGene <- sMetadata %>%
        rstatix::get_summary_stats(nGene, type = "full") %>%
        dplyr::mutate(outlierTreshold_nGene = (q3 + iqr * 1.5)) %>%
        pull(outlierTreshold_nGene) %>%
        unname()

    outlierTreshold_nUMI <- sMetadata %>%
        rstatix::get_summary_stats(nUMI, type = "full") %>%
        dplyr::mutate(outlierTreshold_nUMI = (q3 + iqr * 1.5)) %>%
        pull(outlierTreshold_nUMI) %>%
        unname()

    s_seurat_post_InitFiltering <- subset(
        x = seuratObj,
        subset = (nUMI >= nUMI_l) &
            (nUMI <= outlierTreshold_nUMI) &
            (nGene >= nGene_l) &
            (nGene <= outlierTreshold_nGene) &
            (log10GenesPerUMI > log10GenesPerUMI_l) &
            (mitoRatio < mitoRatio_u)
    )

    s_seurat_post_GeneFiltering <- s_seurat_post_InitFiltering
    counts <- Seurat::GetAssayData(object = s_seurat_post_GeneFiltering, slot = "counts")
    nonzero <- counts > 0
    keepGenes <- Matrix::rowSums(nonzero) >= 10
    filteredCounts <- counts[keepGenes, ]
    return(Seurat::CreateSeuratObject(filteredCounts, meta.data = s_seurat_post_GeneFiltering@meta.data, project = seuratObj@project.name))
}

#' Use publication data to preprocess a Seurat object
#'
#' This function reads in metadata from a CSV file and uses it to preprocess a Seurat object.
#' The function subsets the Seurat object to include only the cells listed in the metadata file,
#' and then performs SCTransform normalization and clustering using the Seurat package.
#'
#' @param seuratObj A Seurat object.
#' @param verbose A logical value indicating whether to print verbose output (default is FALSE).
#' @return A preprocessed Seurat object.
#' @export
UsePublicationData <- function(seuratObj, verbose = FALSE) {
    metadata <- read.csv(file.path("./raw_data", seuratObj@project.name, "metadata.csv"), row.names = 1) |>
        dplyr::filter(cluster_name != "NotAssigned") |>
        dplyr::mutate(UMAP_1 = as.numeric(UMAP_1), UMAP_2 = as.numeric(UMAP_2))
    if ("sample" %in% colnames(seuratObj@meta.data)) {
        metadata <- dplyr::select(metadata, -c(sample))
    }
    seuratObj <- subset(seuratObj, cells = rownames(metadata))
    seuratObj@meta.data <- seuratObj@meta.data |>
        merge(metadata, by = "row.names") |>
        tibble::column_to_rownames("Row.names")
    seuratObj <- Seurat::SCTransform(seuratObj, vst.flavor = "v2", vars.to.regress = "percent.mt", verbose = verbose)
    seuratObj <- GetSeuratClusters(seuratObj, stepName = "usePublicationData")
    Idents(object = seuratObj) <- paste0("SCT_snn_res.", 0.6)
    return(seuratObj)
}
