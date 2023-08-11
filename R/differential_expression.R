#' Get differential expression markers for a given cell type
#'
#' This function takes a Seurat object and returns a list of data frames containing the differential expression markers
#' for a given cell type. The function first subsets the Seurat object by age (if `subsetByAge` is `TRUE`),
#' then loops through each timepoint and cluster to calculate the differential expression markers using the
#' `GetMarkersForCluster` function. The function returns a list of data frames, where each data frame contains the
#' differential expression markers for a single cluster.
#'
#' @param seuratObj A Seurat object.
#' @param cellType A character string specifying the cell type of interest.
#' @param compare A character string specifying the variable to compare between groups (either "experiment" or "sex").
#' @param split A character string specifying the variable to split the comparison by (either "sex", "experiment", or NULL).
#' @param subsetByAge A logical value indicating whether to subset the Seurat object by age.
#' @param logfc.threshold A numeric value specifying the log fold change threshold for differential expression.
#' @param p_val_adj A numeric value specifying the adjusted p-value threshold for differential expression.
#' @param min.pct A numeric value specifying the minimum percentage of cells expressing a gene for it to be considered in the differential expression analysis.
#'
#' @return A list of data frames containing the differential expression markers for each cluster.
#'
#' @export
GetMarkers <- function(seuratObj, cellType, compare = "experiment", split, subsetByAge = FALSE, logfc.threshold = 0.3, p_val_adj = 0.05, min.pct = 0.25) {
    # compare either "Oxygen" or "Male"
    # split either "Sex", "FCG" or NULL
    timepoints <- `if`(subsetByAge, levels(seuratObj$age), age)
    comparisonFolder <- paste0("by", stringr::str_to_title(compare), if (!missing(split)) stringr::str_to_title(split))
    markerList <- list()
    for (time in timepoints) {
        seuratSubset <- subset(seuratObj, subset = age == time)
        for (clusterName in levels(seuratSubset)) {
            markerList[[clusterName]] <- GetMarkersForCluster(seuratSubset,
                cellType = cellType, clusterName = clusterName,
                split = split, compare = compare, time = time, logfc.threshold = logfc.threshold, FDR = p_val_adj, min.pct = min.pct
            )
        }
    }
    return(markerList)
}

#' Get differential expression markers for a given cluster
#'
#' This function takes a Seurat object subset and returns a data frame containing the differential expression markers for a given cluster.
#' The function first subsets the Seurat object by cluster and (optionally) by a given split variable. The function then performs differential
#' expression analysis using the `FindMarkers` function and returns a data frame containing the differential expression markers for the given cluster.
#'
#' @param seuratSubset A Seurat object.
#' @param cellType A character string specifying the cell type of interest.
#' @param clusterName A character string specifying the name of the cluster of interest.
#' @param split A character string specifying the variable to split the comparison by (either "sex", "experiment", or NULL).
#' @param compare A character string specifying the variable to compare between groups (either "experiment" or "sex").
#' @param time A character string specifying the timepoint of the Seurat object subset.
#' @param logfc.threshold A numeric value specifying the log fold change threshold for differential expression.
#' @param FDR A numeric value specifying the adjusted p-value threshold for differential expression.
#' @param min.pct A numeric value specifying the minimum percentage of cells expressing a gene for it to be considered in the differential expression analysis.
#'
#' @return A data frame containing the differential expression markers for the given cluster.
#'
#' @export
GetMarkersForCluster <- function(seuratSubset, cellType, clusterName, split, compare, time, logfc.threshold, FDR, min.pct) {
    comparisonFolder <- paste0("by", stringr::str_to_title(compare), if (!missing(split)) stringr::str_to_title(split))
    splitBy <- c("None")
    if (!missing(split)) {
        splitBy <- levels(seuratSubset@meta.data[[split]])
    }
    xyDf <- data.frame()
    DEresults <- list()
    for (type in splitBy) {
        tryCatch(
            {
                seu <- subset(x = seuratSubset, idents = ifelse(clusterName == "acap_reactive_acap", c("aCap", "Reactive aCap"), clusterName))
                if (type != "None") {
                    seu <- seu[, seu@meta.data[, split] == type]
                }
                markers <- Seurat::SetIdent(seu, value = compare) |>
                    Seurat::FindMarkers(
                        assay = "SCT", ident.1 = levels(seuratSubset@meta.data[[compare]])[2],
                        ident.2 = levels(seuratSubset@meta.data[[compare]])[1],
                        logfc.threshold = 0.2, recorrect_umi = FALSE
                    ) |>
                    dplyr::mutate(sig_exp = dplyr::case_when(
                        avg_log2FC > logfc.threshold & p_val_adj < FDR ~ "UP",
                        avg_log2FC < -logfc.threshold & p_val_adj < FDR ~ "DOWN",
                        .default = "NO"
                    ))

                xyCh <- openxlsx::read.xlsx(file.path("data", "XY_ch.xlsx"))

                markers <- markers %>%
                    tibble::rownames_to_column(var = "Gene") %>%
                    dplyr::left_join(xyCh, by = "Gene") %>%
                    dplyr::relocate("Chromosome", .after = "Gene") %>%
                    tibble::column_to_rownames(var = "Gene")

                if (type == "None") {
                    DEresults <- markers
                } else {
                    DEresults[[type]] <- markers
                }

                dir.create(file.path("./results", age, "byCellType", cellType, comparisonFolder, clusterName),
                    recursive = TRUE, showWarnings = FALSE
                )

                openxlsx::write.xlsx(markers,
                    file = file.path(
                        "./results", age, "byCellType", cellType, comparisonFolder,
                        clusterName, paste0(
                            clusterName, if (type != "None") paste0("_", type),
                            "_", levels(seuratSubset@meta.data[[compare]])[2], "_vs_",
                            levels(seuratSubset@meta.data[[compare]])[1], ".xlsx"
                        )
                    ), quote = FALSE, rowNames = TRUE
                )

                xyGenes <- markers |>
                    dplyr::filter(sig_exp != "NO" & !is.na(Chromosome)) |>
                    dplyr::select(Chromosome, sig_exp) |>
                    tibble::rownames_to_column(var = "Gene") |>
                    dplyr::mutate(cluster = clusterName)
                if (type != "None") {
                    xyGenes <- xyGenes |> dplyr::mutate(type = type)
                }
                xyDf <- rbind(xyDf, xyGenes)
            },
            error = function(e) {
                cat("ERROR :", conditionMessage(e), "\n")
            }
        )
    }
    if ("type" %in% colnames(xyDf)) {
        xyDf <- xyDf |>
            dplyr::group_by(cluster, type, Chromosome, sig_exp)
    } else {
        xyDf <- xyDf |>
            dplyr::group_by(cluster, Chromosome, sig_exp)
    }
    xyDf <- xyDf |> dplyr::summarise(number = dplyr::n())
    write.csv(xyDf, file = file.path("./results", time, "byCellType", cellType, comparisonFolder, "XY_number.csv"))
    return(DEresults)
}
