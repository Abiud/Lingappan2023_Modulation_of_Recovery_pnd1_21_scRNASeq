#' Get Seurat markers
#'
#' This function takes a Seurat object and returns a list of data frames containing the top markers for each cluster,
#' as well as a spreadsheet containing additional information about the markers.
#'
#' @param seuratObj A Seurat object.
#'
#' @return A list of data frames containing the top markers for each cluster.
#' @examples
#' seuratObj <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19/")
#' seuratObj <- CreateSeuratObject(counts = seuratObj)
#' seuratObj <- NormalizeData(seuratObj)
#' seuratObj <- FindVariableFeatures(seuratObj)
#' seuratObj <- ScaleData(seuratObj)
#' seuratObj <- RunPCA(seuratObj)
#' seuratObj <- FindNeighbors(seuratObj)
#' seuratObj <- FindClusters(seuratObj)
#' markers <- GetSeuratMarkers(seuratObj)
#' @export
GetSeuratMarkers <- function(seuratObj) {
    # set future options
    options(future.globals.maxSize = 31457280000)
    future::plan("multicore", workers = future::availableCores() - 1)
    outputName <- paste0(age, "_markers")
    refSheets <- c("pnd1", "pnd21")
    refList <- list()
    for (refName in refSheets) {
        refList[refName] <- openxlsx::read.xlsx(paste0("./data/cluster_annotation_files/", refName, "_clusters.xlsx")) |>
            list()
    }

    # styling
    headerDecoration <- openxlsx::createStyle(
        textDecoration = "Bold",
        fontSize = 14, fontColour = "#FFFFFF", halign = "center",
        fgFill = "#4F81BD", border = "TopBottom", borderColour = "#4F81BD"
    )
    not0Style <- openxlsx::createStyle(fontColour = "#000000", bgFill = "#EDFFF0")
    g5Style <- openxlsx::createStyle(fontColour = "#000000", bgFill = "#D0F7F7")
    g10Style <- openxlsx::createStyle(fontColour = "#000000", bgFill = "#FAC3C3")
    # start an empty dataframe
    dfHolder <- data.frame(
        Date = as.Date(character()),
        File = character(),
        User = character(),
        stringsAsFactors = FALSE
    )

    # get empty data.frame that will store all find markers
    dfAll <- dfHolder
    totalRows <- 0

    refCountList <- list()
    for (refName in refSheets) {
        refCount <- data.frame(matrix(ncol = ncol(refList[[refName]]) + 2, nrow = nlevels(seuratObj)))
        colnames(refCount) <- c("cluster_name", colnames(refList[[refName]]), "possible_name")
        refCount$cluster_name <- levels(seuratObj)
        refCountList[refName] <- list(refCount)
    }

    clusterMarkerList <- list()
    wb <- openxlsx::createWorkbook()
    for (clusterName in levels(seuratObj)) {
        foundMarkers <- Seurat::FindMarkers(seuratObj,
            ident.1 = clusterName, verbose = FALSE,
            only.pos = TRUE, test.use = "wilcox", assay = "SCT"
        )
        # sort markers by avg_log2FC
        df <- dplyr::arrange(foundMarkers, desc(foundMarkers$avg_log2FC))
        clusterMarkerList[[clusterName]] <- df
        # Save cluster worksheet
        sheetName <- paste0("Cluster_", clusterName)
        if (nchar(sheetName) > 31) {
            sheetName <- paste0("Cluster_", shortenClusterName(clusterName))
        }
        openxlsx::addWorksheet(wb,
            sheetName = sheetName,
            gridLines = TRUE, zoom = 180
        )
        openxlsx::writeData(wb,
            sheet = sheetName, x = df,
            headerStyle = headerDecoration, rowNames = TRUE
        )
        # Add markers to list of all markers
        df$cluster_name <- clusterName
        df <- df |>
            tibble::rownames_to_column(var = "gene") |>
            dplyr::relocate(gene, .before = p_val)
        row.names(df) <- totalRows + seq_len(nrow(df))
        totalRows <- totalRows + nrow(df)
        dfAll <- rbind.data.frame(dfAll, df)

        # get top 25 markers of current cluster and loop trough columns of the already classified markers
        # to check how many intersect on each, if they match add avg_log2FC score to score table
        tmpDf <- head(df, 25)
        rownames(tmpDf) <- tmpDf$gene
        for (refName in refSheets) {
            for (clusterN in colnames(refList[[refName]])) {
                tmpInt <- dplyr::intersect(refList[[refName]][, clusterN], tmpDf[, "gene"])
                avg2Sum <- 0
                for (item in tmpInt) {
                    avg2Sum <- avg2Sum + tmpDf[item, "avg_log2FC"]
                }
                refCountList[[refName]][refCountList[[refName]]$cluster_name == clusterName, clusterN] <- round(avg2Sum, 2)
            }
            refCountList[[refName]][refCountList[[refName]]$cluster_name == clusterName, "possible_name"] <-
                refCountList[[refName]][refCountList[[refName]]$cluster_name == clusterName, ] %>%
                tidyr::gather(col, val, dplyr::first(colnames(refList[[refName]])):dplyr::last(colnames(refList[[refName]]))) %>%
                dplyr::filter(val == max(val)) %>%
                dplyr::select(-cluster_name, -val, -possible_name)
        }
    }

    openxlsx::addWorksheet(wb, sheetName = "all_clusters", gridLines = TRUE, zoom = 180)
    openxlsx::freezePane(wb, sheet = "all_clusters", firstRow = TRUE, firstCol = TRUE)
    openxlsx::writeData(wb,
        sheet = "all_clusters", x = dfAll,
        headerStyle = headerDecoration, rowNames = TRUE
    )

    openxlsx::addWorksheet(wb, sheetName = "top100_common", gridLines = TRUE, zoom = 180)
    openxlsx::freezePane(wb, sheet = "top100_common", firstCol = TRUE)

    #####
    startRow <- 1
    for (refName in refSheets) {
        openxlsx::mergeCells(wb, "top100_common", cols = 2:8, rows = startRow)
        openxlsx::writeData(wb, "top100_common",
            paste0(
                "Comparing top 25 genes on ",
                refName,
                " to found top 25 genes on each cluster"
            ),
            startCol = 2, startRow = startRow,
        )
        openxlsx::conditionalFormatting(wb, "top100_common",
            cols = 2:ncol(refCountList[[refName]]),
            rows = (startRow + 2):(nrow(refCountList[[refName]]) + startRow + 1),
            rule = "!=0", style = not0Style
        )
        openxlsx::conditionalFormatting(wb, "top100_common",
            cols = 2:ncol(refCountList[[refName]]),
            rows = (startRow + 2):(nrow(refCountList[[refName]]) + startRow + 1),
            rule = ">=5", style = g5Style
        )
        openxlsx::conditionalFormatting(wb, "top100_common",
            cols = 2:ncol(refCountList[[refName]]),
            rows = (startRow + 2):(nrow(refCountList[[refName]]) + startRow + 1),
            rule = ">=10", style = g10Style
        )
        openxlsx::addStyle(wb, "top100_common", cols = seq_len(ncol(refCountList[[refName]])), rows = startRow + 1, style = headerDecoration)
        openxlsx::writeData(wb,
            sheet = "top100_common", x = refCountList[[refName]],
            headerStyle = headerDecoration, rowNames = FALSE, startRow = startRow + 1
        )

        startRow <- startRow + nrow(refCountList[[refName]]) + 3
    }

    # Get unique markers on each cluster
    # create a list containing the lists of all clusters marker names
    clusterList <- list()
    for (clusterName in levels(seuratObj)) {
        sheetName <- paste0("Cluster_", clusterName)
        if (nchar(sheetName) > 31) {
            sheetName <- paste0("Cluster_", shortenClusterName(clusterName))
        }
        clusterList[[clusterName]] <- rownames(openxlsx::read.xlsx(wb, sheet = sheetName, rowNames = TRUE))
    }

    # compare every cluster against all other cluster to discover unique markers
    for (i_clusterName in levels(seuratObj)) {
        df <- head(clusterList[[i_clusterName]], 100)

        for (j_clusterName in levels(seuratObj)) {
            if (i_clusterName == j_clusterName) next
            equal <- dplyr::intersect(df, clusterList[[j_clusterName]])
            df <- dplyr::setdiff(df, equal)
        }
        sheetName <- paste0("Cluster_", i_clusterName)
        if (nchar(sheetName) > 31) {
            sheetName <- paste0("Cluster_", shortenClusterName(i_clusterName))
        }
        result <- openxlsx::read.xlsx(wb, sheet = sheetName, rowNames = TRUE)
        result <- result[df, ]

        sheetName <- paste0("unique_cl_", i_clusterName)
        if (nchar(sheetName) > 31) {
            sheetName <- paste0("unique_cl_", shortenClusterName(i_clusterName))
        }
        openxlsx::addWorksheet(wb,
            sheetName = sheetName,
            gridLines = TRUE, zoom = 180
        )
        openxlsx::writeData(wb,
            sheet = sheetName, x = result,
            headerStyle = headerDecoration, rowNames = TRUE
        )
    }
    dir.create(file.path("results", age), recursive = TRUE)
    openxlsx::saveWorkbook(wb,
        file = paste0("./results/", age, "/", outputName, ".xlsx"),
        overwrite = TRUE
    )

    return(clusterMarkerList)
}
