#' Integrate multiple Seurat objects
#'
#' This function integrates multiple seurat objects
#'
#' @param seuratList List of Seurat objects to be integrated
#' @param fiaKAnchor k.anchor parameter on Seurat::FindIntegrationAnchors, default : 5
#' @param fiaReduction reduction method on Seurat::FindIntegrationAnchors, default : "cca"
#' @return Integrated Seurat Object
#' @export
AggregateSamples <- function(seuratList, fiaKAnchor = 5, fiaReduction = "cca", verbose = FALSE) {
    # set future options
    options(future.globals.maxSize = 31457280000)
    future::plan("multicore", workers = future::availableCores() - 1)
    # Start integration Pipeline using the list of Seurat objects
    seuratFeatures <- Seurat::SelectIntegrationFeatures(object.list = seuratList, nfeatures = 3000, verbose = verbose)
    seuratList <- Seurat::PrepSCTIntegration(
        object.list = seuratList, anchor.features = seuratFeatures,
        verbose = verbose
    )
    # Next, identify anchors and integrate the datasets.
    # make sure to set normalization.method = 'SCT'
    seuratAnchor <- Seurat::FindIntegrationAnchors(
        object.list = seuratList, normalization.method = "SCT",
        anchor.features = seuratFeatures, verbose = verbose, reference = 1,
        k.anchor = fiaKAnchor, reduction = fiaReduction
    )
    seuratIntegrated <- Seurat::IntegrateData(
        anchorset = seuratAnchor, normalization.method = "SCT",
        verbose = verbose
    )

    seuratIntegrated@project.name <- paste0("integrated_", age)
    dir.create(file.path("figures", age, "QC", seuratIntegrated@project.name), recursive = TRUE, showWarnings = FALSE)
    seuratIntegrated <- GetSeuratClusters(seuratIntegrated, TSNE = FALSE, stepName = "IntegrateData")

    seuratIntegrated@meta.data <- seuratIntegrated@meta.data %>%
        dplyr::mutate(
            sex = dplyr::case_when(
                grepl("F", sample) ~ "Female",
                grepl("M", sample) ~ "Male",
            ),
            experiment = dplyr::case_when(
                grepl("RoomAir", sample) ~ "Room_Air",
                grepl("Oxygen", sample) ~ "Oxygen"
            ),
            age = age
        )

    seuratIntegrated$sex <- factor(seuratIntegrated$sex, levels = c("Female", "Male"))
    seuratIntegrated$experiment <- factor(seuratIntegrated$experiment, levels = c("Room_Air", "Oxygen"))
    seuratIntegrated$celltype <- factor(seuratIntegrated$celltype)
    seuratIntegrated$cluster_name <- factor(seuratIntegrated$cluster_name)
    seuratIntegrated$age <- factor(seuratIntegrated$age)
    seuratIntegrated$sample <- factor(seuratIntegrated$sample)
    seuratIntegrated$library_name <- factor(seuratIntegrated$library_name)

    umap_mat <- seuratIntegrated@meta.data |>
        dplyr::select(c(UMAP_1, UMAP_2)) |>
        as.matrix()
    seuratIntegrated[["UMAP"]] <- Seurat::CreateDimReducObject(umap_mat, key = "UMAP_", assay = "SCT", global = TRUE)

    Seurat::DefaultAssay(seuratIntegrated) <- "SCT"
    seuratIntegrated <- Seurat::PrepSCTFindMarkers(seuratIntegrated, verbose = verbose)
    Idents(seuratIntegrated) <- "cluster_name"
    return(seuratIntegrated)
}
