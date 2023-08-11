#' Set the SoupX channel for a sample
#'
#' @param sampleDir A character string specifying the directory containing the sample data
#'
#' @return A data frame containing the estimated SoupX channel
#'
#' @export
SetSoupChannel <- function(sampleDir) {
    toc <- Seurat::Read10X(file.path(sampleDir))
    tod <- Seurat::Read10X(file.path(sampleDir, "raw_feature_bc_matrix"))
    sc <- SoupX::SoupChannel(tod, toc, calcSoupProfile = FALSE)
    return(SoupX::estimateSoup(sc))
}

#' Get adjusted counts for a Seurat sample using SoupX
#'
#' @param seuratSample A Seurat object containing gene expression data for a specific sample
#' @param soupChannel A data frame containing the raw counts for the SoupX channel
#' @param contamination A numeric value specifying the estimated contamination fraction, or "auto" to estimate automatically
#'
#' @return A list containing the adjusted counts and the estimated contamination fraction
#'
#' @export
GetSoupAdjustedCounts <- function(seuratSample, soupChannel) {
    sc_auto <- SoupX::setContaminationFraction(soupChannel, 0.05)
    sc <- SoupX::setClusters(sc_auto, setNames(seuratSample$SCT_snn_res.0.6, rownames(sc_auto)))
    return(SoupX::adjustCounts(sc))
}

#' Set the SoupX channel
#'
#' This function sets the SoupX channel for a given sample.
#'
#' @param sample A Seurat object.
#'
#' @return  A data frame containing the estimated SoupX channel.
#'
#' @export
GetSoupChannel <- function(sample) {
    return(SetSoupChannel(paste0("./raw_data/", sample@project.name)))
}

#' Run SoupX on a sample
#'
#' This function runs SoupX on a given sample and returns a preprocessed Seurat object.
#'
#' @param sample A Seurat object.
#' @param soupChannel A character string representing the SoupX channel.
#'
#' @return A preprocessed Seurat object.
#'
#' @export
RunSoupX <- function(sample, soupChannel) {
    sampleName <- sample@project.name
    soupResults <- GetSoupAdjustedCounts(sample, soupChannel)
    postSoupSample <- CreateSeuratObject(counts = soupResults, project = sampleName)
    postSoupSample$percent.mt <- Seurat::PercentageFeatureSet(object = postSoupSample, pattern = "^mt-")
    postSoupSample$percent.ribo <- Seurat::PercentageFeatureSet(postSoupSample, pattern = "^Rp[sl]")
    postSoupSample@meta.data <- postSoupSample@meta.data |>
        dplyr::mutate(
            log10GenesPerUMI = log10(nFeature_RNA) / log10(nCount_RNA),
            mitoRatio = percent.mt / 100,
            riboRatio = percent.ribo / 100,
            nGene = nFeature_RNA,
            nUMI = nCount_RNA,
            sample = orig.ident
        )

    return(postSoupSample)
}
