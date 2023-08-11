#' Detect doublets in a Seurat object
#'
#' This function detects doublets in a given Seurat object using the DoubletFinder package.
#' The function first calculates the optimal parameters for doublet detection using the
#' `GetDoubletParams` function, then performs doublet detection using the `doubletFinder_v3`
#' function. The function also performs a second round of doublet detection using a lower doublet
#' formation rate to remove any remaining doublets.
#'
#' @param seuratObj A Seurat object.
#'
#' @return A preprocessed Seurat object with doublets removed.
#'
#' @export
GetDoublets <- function(seuratObj, dParams) {
    seuratObj <- DoubletFinder::doubletFinder_v3(seuratObj,
        PCs = 1:dParams$dims, pK = dParams$optimalPK, nExp = dParams$nExp_poi,
        reuse.pANN = FALSE, sct = TRUE
    )

    pANN <- paste("pANN", 0.25, dParams$optimalPK, dParams$nExp_poi, sep = "_")
    seuratObj <- DoubletFinder::doubletFinder_v3(seuratObj,
        PCs = 1:30, pK = dParams$optimalPK, nExp = dParams$nExp_poi.adj,
        reuse.pANN = pANN,
        sct = TRUE
    )
    return(seuratObj)
}

#' Get doublet detection parameters for a Seurat object
#'
#' This function calculates the optimal parameters for doublet detection using the DoubletFinder package.
#' The function first calculates the number of significant principal components to use for doublet detection,
#' then performs a parameter sweep to find the optimal bimodality coefficient threshold. Finally,
#' the function calculates the expected number of doublets based on the number of cells in the sample and the homotypic doublet rate.
#'
#' @param seuratObj A Seurat object.
#'
#' @return A list containing the following elements:
#' - dims: the number of significant principal components to use for doublet detection
#' - optimalPK: the optimal bimodality coefficient threshold
#' - nExp_poi: the expected number of doublets based on the number of cells in the sample and a 7.5% doublet formation rate
#' - nExp_poi.adj: the expected number of doublets adjusted for homotypic doublet rate
#'
#' @export
GetDoubletParams <- function(seuratObj) {
    dimsCalc <- GetSignificantPC(seuratObj, "DoubletFindAndRemove", age)
    invisible(capture.output(sweepResList <- DoubletFinder::paramSweep_v3(seuratObj, PCs = 1:dimsCalc, sct = TRUE, num.cores = future::availableCores())))
    sweepStats <- DoubletFinder::summarizeSweep(sweepResList)
    invisible(capture.output(bcmvn_s_sct <- DoubletFinder::find.pK(sweepStats))) # bcmvn -> mean-variance normalized bimodality coefficient
    optimalPK <- as.numeric(as.character(bcmvn_s_sct$pK[which.max(bcmvn_s_sct$BCmetric)]))

    annotations <- seuratObj@meta.data$seurat_clusters # chosen.resolution = 0.6
    homotypic.prop <- DoubletFinder::modelHomotypic(annotations) ## ex: annotations <- seuratObj@meta.data$seurat_clusters
    nExp_poi <- round(0.075 * nrow(seuratObj@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

    return(list(dims = dimsCalc, optimalPK = optimalPK, nExp_poi = nExp_poi, nExp_poi.adj = nExp_poi.adj))
}

#' Remove doublets from a Seurat object
#'
#' This function removes doublets from a given Seurat object using the DoubletFinder package.
#' The function first selects the singlet cells using the `DF.classifications` metadata column,
#' then performs clustering using the `GetSeuratClusters` function. Finally, the function sets the
#' cluster identities using the `Idents` function and returns the preprocessed Seurat object with doublets removed.
#'
#' @param seuratObj A Seurat object.
#' @param dParams A list containing the doublet detection parameters.
#'
#' @return A preprocessed Seurat object with doublets removed.
#'
#' @export
#'
#' @examples
#' seuratObj <- Read10X(data.dir = "data/filtered_gene_bc_matrices/")
#' seuratObj <- CreateSeuratObject(counts = seuratObj)
#' doubletParams <- GetDoubletParams(seuratObj)
#' seuratObj <- GetDoublets(seuratObj, doubletParams)
#' seuratObj <- RemoveDoublets(seuratObj, doubletParams)
RemoveDoublets <- function(seuratObj, dParams) {
    cells.use <- colnames(seuratObj)[with(
        seuratObj@meta.data,
        get(paste("DF.classifications", 0.25, dParams$optimalPK, dParams$nExp_poi.adj, sep = "_"))
    ) == "Singlet"]

    seuratObj <- GetSeuratClusters(subset(seuratObj, cells = cells.use), stepName = "removeCont2")
    Idents(object = seuratObj) <- paste0("SCT_snn_res.", 0.6) # "SCT_snn_res.0.6"
    return(seuratObj)
}
