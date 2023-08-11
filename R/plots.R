GetIntegratedPlots <- function(seuratObj) {
    dimplots <- GetDimPlots(seuratObj)
    featurePlots <- Seurat::FeaturePlot(seuratObj, features = c("Epcam", "Cldn5", "Ptprc", "Col1a2"))

    dotPlot <- GenerateDotPlot(seuratObj, c(
        "Gpihbp1", "Car4", "Cxcl12", "Vwf", "Top2a", "Ager", "Cldn18",
        "Retnla", "Gal", "C1qa", "H2-Aa", "Spp1", "Ly6c2", "Stfa1",
        "Nr4a1", "Retnlg", "Ngp", "Ighm", "Gzma", "Gucy1a1", "Cd200r3",
        "Tagln", "Mfap4"
    ))

    return(list(dimPlots = dimplots, featurePlots = featurePlots, dotPlot = dotPlot)) # , featurePlots = featurePlots, dotPlot = dotPlot))
}


GetSubsetPlots <- function(seuratObj) {
    dimPlots <- GetDimPlots(seuratObj)
    proportionPlots <- GetProportionPlots(seuratObj)
    return(list(dimPlots = dimPlots, proportionPlots = proportionPlots))
}

#' Generate and save dimensionality reduction plots for a Seurat object
#'
#' @param seuratObj A Seurat object containing gene expression data for a specific sample
#' @param cellType A character string specifying the cell type to plot
#' @param colorPalette A character vector specifying the color palette to use for the plots
#' @param FCG A logical value specifying whether to plot by sex_ch and sex_gon
#'
#' @return A list of ggplot objects containing the generated plots
#'
#' @export
GetDimPlots <- function(seuratObj, colorPalette = c(
                            "#F8766D", "#3DA1FF", "#DD8D00", "#BE80FF", "#B3A000",
                            "#97A900", "#FF6C92", "#2FB600", "#8F91FF", "#00BF76",
                            "#FE61CF", "#00C0B7", "#00BDD1", "#F265E7", "#00AEFA",
                            "#EC823C", "#00BB4B", "#CA9700", "#DE71F9", "#00B7E8",
                            "#00C098", "#FF64B3", "#71B000", "#FF64B0"
                        )) {
    plots <- list()
    umap <- Seurat::DimPlot(seuratObj, label = TRUE, cols = colorPalette[seq_len(nlevels(seuratObj))])
    splitBySample <- Seurat::DimPlot(seuratObj, split.by = "sample", ncol = 2, label = TRUE, cols = colorPalette[seq_len(nlevels(seuratObj))]) +
        Seurat::NoLegend()
    splitByExperiment <- Seurat::DimPlot(seuratObj, split.by = "experiment", cols = colorPalette[seq_len(nlevels(seuratObj))])
    groupBySex <- Seurat::DimPlot(seuratObj, group.by = "sex", split.by = "experiment", cols = c("#9460A0", "#97CED0"), shuffle = TRUE)

    plots[["umap"]] <- umap
    plots[["splitBySample"]] <- splitBySample
    plots[["splitByExperiment"]] <- splitByExperiment
    plots[["groupBySex"]] <- groupBySex

    return(plots)
}

#' Generate and save proportion plots for a Seurat object
#'
#' @param seuratObj A Seurat object containing gene expression data for a specific sample
#' @param cellType A character string specifying the cell type to plot
#' @param colorPalette A character vector specifying the color palette to use for the plots
#'
#' @return A list of ggplot objects containing the generated plots
#'
#' @export
GetProportionPlots <- function(seuratObj, colorPalette) {
    if (missing(colorPalette)) {
        colorPalette <- ggplotColours(length(levels(seuratObj)))
    }
    yLevels <- levels(seuratObj)

    pt <- table(Seurat::Idents(seuratObj), seuratObj$experiment)
    pt <- as.data.frame(pt)
    pt$Var1 <- as.character(pt$Var1)
    pt$Var2 <- factor(pt$Var2, levels = c("Room_Air", "Oxygen"))
    plotA <- ggplot2::ggplot(pt, ggplot2::aes(x = Var2, y = Freq, fill = factor(Var1, levels = yLevels))) +
        ggplot2::theme_bw(base_size = 15) +
        ggplot2::geom_col(position = "fill", width = 0.5) +
        ggplot2::xlab("Sample") +
        ggplot2::ylab("Proportion") +
        ggplot2::scale_y_continuous(labels = scales::percent) +
        ggplot2::scale_fill_manual(values = colorPalette) +
        ggplot2::theme_bw()

    pt <- table(Seurat::Idents(seuratObj), seuratObj$sample)
    pt <- as.data.frame(pt)
    pt$Var1 <- as.character(pt$Var1)
    plotB <- ggplot2::ggplot(pt, ggplot2::aes(x = Var2, y = Freq, fill = factor(Var1, levels = yLevels))) +
        ggplot2::geom_col(position = "fill", width = 0.5) +
        ggplot2::xlab("Sample") +
        ggplot2::scale_fill_manual(values = colorPalette) +
        ggplot2::scale_y_continuous(labels = scales::percent) +
        ggplot2::theme_bw()

    return(list(bySample = plotA, byExperiment = plotB))
}

#' Generate Venn diagram
#'
#' This function generates the venn diagrams for the given data
#'
#' @param data named list with two vectors to plot
#' @param clusterName cell type of object being analyzed
#' @param labels wheter venn diagram should include labels or not
#' @param colors vector with two hex colors to use
#' @return ggplot2 object
#' @export
GetVenn <- function(data, clusterName, labels = FALSE, colors = c("#97CED0", "#9460A0")) {
    return(VennDiagram::venn.diagram(
        x = list(data[[1]], data[[2]]), category.names = names(data),
        filename = NULL, disable.logging = TRUE,
        rotation.degree = ifelse(length(data[[1]]) < length(data[[2]]), 180, 0),
        scaled = FALSE, main = clusterName, main.fontface = "bold", main.fontfamily = "sans",
        # Circles
        lwd = 2,
        lty = "blank",
        fill = colors,
        # Numbers
        fontface = "bold", fontfamily = "sans",
        # Set names
        cat.fontface = "bold", cat.default.pos = "outer", cat.fon = "sans",
        width = 8, height = 8, units = "in"
    ))
}

#' Generate Heatmaps
#'
#' This function generates the venn diagrams for each seurat cluster and
#' saves an excel report about the unique and shared genes.
#'
#' @param seuratObj Seurat object being analyzed
#' @param cellType cell type of object being analyzed
#' @param experiment type of comparison
#' @param suffix suffix of clusters marker files
#' @return Null
#' @export
GenerateHeatmaps <- function(seuObj, clusterName, DEmarkers, compareBy = "sex", experiment = "bySex") {
    compareFactors <- switch(compareBy,
        "sex" = c("Male", "Female"),
        "experiment" = c("Oxygen", "Room_Air"),
        "experimentSex" = c("Oxygen", "Room_Air"),
        NA
    )
    seu <- subset(seuObj, idents = clusterName)
    if (ncol(seu) < 50) {
        print(paste("Skipping", clusterName, "because it has less than 50 cells"))
        return()
    }
    if (compareBy == "experimentSex") {
        if (!is.null(DEmarkers$Male)) {
            getHeatmaps(
                seu, DEmarkers$Male, "experiment", experiment, compareFactors, clusterName, "Male"
            )
        }
        if (!is.null(DEmarkers$Female)) {
            getHeatmaps(
                seu, DEmarkers$Female, "experiment", experiment, compareFactors, clusterName, "Female"
            )
        }
    } else {
        getHeatmaps(seu, DEmarkers, compareBy, experiment, compareFactors, clusterName)
    }
}

getHeatmaps <- function(seuObj, DEmarkers, compareBy, experiment, compareFactors, clusterName, sexName) {
    if (!missing(sexName)) {
        seuObj <- subset(seuObj, subset = sex == sexName)
    }
    sexColorPalette <- c("#97CED0", "#9460A0")
    experimentColorPalette <- c("#009292", "#E2B038")
    compCol <- unlist(seuObj[[compareBy]])
    names(compCol) <- rownames(seuObj[[compareBy]])
    compCol <- factor(compCol, compareFactors)
    seuObj[[compareBy]] <- compCol
    seuObj <- Seurat::SetIdent(seuObj, value = compareBy)
    seuObj <- Seurat::SCTransform(seuObj, vst.flavor = "v2", verbose = FALSE, return.only.var.genes = FALSE)
    seuObj <- Seurat::RunPCA(seuObj, verbose = FALSE)

    df <- DEmarkers |>
        dplyr::arrange(desc(avg_log2FC)) |>
        dplyr::filter(sig_exp %in% c("UP", "DOWN"))

    if (nrow(df) < 5) {
        print(paste(clusterName, "skipped because it has less than 5 DE genes"))
        return()
    }

    seuMeta <- seuObj@meta.data |>
        dplyr::select(sex) |>
        dplyr::arrange(sex)
    if (compareBy == "experiment") {
        seuMeta <- seuObj@meta.data |>
            dplyr::select(experiment, sex) |>
            dplyr::arrange(experiment, desc(sex))
        oxyCells <- rownames(seuMeta[seuMeta$experiment == "Oxygen", ])
        roomCells <- rownames(seuMeta[seuMeta$experiment == "Room_Air", ])
    }
    femaleCells <- rownames(seuMeta[seuMeta$sex == "Female", , drop = FALSE])
    maleCells <- rownames(seuMeta[seuMeta$sex == "Male", , drop = FALSE])

    scaleData <- seuObj@assays$SCT@scale.data[rownames(seuObj@assays$SCT@scale.data) %in% rownames(df), , drop = FALSE] |>
        as.data.frame() |>
        tibble::rownames_to_column(var = "gene") |>
        dplyr::arrange(match(gene, rownames(df))) |>
        tibble::column_to_rownames(var = "gene") |>
        dplyr::select(dplyr::all_of(rownames(seuMeta)))

    # Col Annotation
    colAnnList <- list()
    colColors <- list()
    if (compareBy == "experiment") {
        factorList <- list()
        colAnnList[[compareBy]] <- seuMeta$experiment
        factorList[[compareFactors[1]]] <- experimentColorPalette[1]
        factorList[[compareFactors[2]]] <- experimentColorPalette[2]
        colColors[[compareBy]] <- unlist(factorList)
    }
    factorList <- list()
    colAnnList[["sex"]] <- seuMeta$sex
    factorList[["Male"]] <- sexColorPalette[1]
    factorList[["Female"]] <- sexColorPalette[2]
    colColors[["sex"]] <- unlist(factorList)

    colAnnList[["col"]] <- colColors
    col_ann <- do.call(ComplexHeatmap::HeatmapAnnotation, colAnnList)
    # Row Annotation - select rownames to display
    topNGenes <- df |>
        dplyr::filter(sig_exp == "UP") |>
        dplyr::slice_head(n = 5) |>
        rownames()
    botNGenes <- df |>
        dplyr::filter(sig_exp == "DOWN") |>
        dplyr::slice_tail(n = 5) |>
        rownames()
    rowsToAnnotate <- c(topNGenes, botNGenes)
    ha <- ComplexHeatmap::rowAnnotation(foo = ComplexHeatmap::anno_mark(
        at = which(rownames(scaleData) %in% rowsToAnnotate),
        labels = rowsToAnnotate
    ))
    # Heatmap
    ht <- ComplexHeatmap::Heatmap(scaleData,
        border = TRUE,
        cluster_rows = FALSE,
        show_column_names = FALSE,
        cluster_columns = FALSE,
        name = "Expression",
        column_title = paste0("Cells (", ncol(scaleData), ")"),
        row_title = "DE genes",
        top_annotation = col_ann,
        show_row_names = FALSE,
        right_annotation = ha,
        col = circlize::colorRamp2(
            c(-1, -0.5, 0, 0.5, 1),
            rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu"))
        ),
        use_raster = TRUE
    )

    ComplexHeatmap::draw(ht)
}

#' Generate DotPlot
#'
#' This function generates the DotPlot of the levels in the given seurat object
#'
#' @param seuratObj Seurat object being analyzed
#' @param markersToPlot list of markers to plot
#' @return Null
#' @export
GenerateDotPlot <- function(seuratObj, markersToPlot, dot.scale = 8) {
    levels(seuratObj) <- rev(levels(seuratObj))
    dotPlot <- Seurat::DotPlot(seuratObj, features = markersToPlot, cols = c("#28A4B3", "#C51517"), dot.scale = dot.scale) +
        RotatedAxis()
    return(dotPlot)
}
