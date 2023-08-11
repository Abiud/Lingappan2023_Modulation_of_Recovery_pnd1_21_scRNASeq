#' Run Monocle analysis
#'
#' This function performs Monocle analysis on a Seurat object and generates various plots and analysis results related to the cell trajectory and gene modules.
#'
#' @param seuratObj A Seurat object containing single-cell RNA-seq data.
#' @param age The age category to consider for analysis.
#' @param cellType The cell type to consider for analysis.
#' @param resolution The resolution values to use for community detection during analysis. Defaults to a sequence of powers of 10 ranging from 10^-6 to 10^-1.
#' @param assay The assay to use for analysis. Defaults to "SCT".
#' @param useMonocleCluster Logical value indicating whether to use Monocle cluster assignments or Seurat cluster assignments during analysis. Defaults to TRUE.
#'
#' @return None.
#'
#' @examples
#' RunMonocle(seuratObj, age = "pnd7", cellType = "Alveolar Macrophages", resolution = 10^seq(-6, -1), assay = "SCT", useMonocleCluster = TRUE)
#'
#' @export
RunMonocle <- function(seuratObj, cellType, resolution = 10^seq(-6, -1), assay = "SCT", useMonocleCluster = TRUE) {
    # other resolution <- 1e-03
    saveDir <- paste0("./data/monocle/", cellType)
    figDir <- file.path(saveDir, "figures")
    dataDir <- file.path(saveDir, "data")
    dir.create(figDir, showWarnings = FALSE, recursive = TRUE)
    dir.create(dataDir, showWarnings = FALSE)

    if (!file.exists(file.path(saveDir, "cds.RData"))) {
        DefaultAssay(seuratObj) <- assay
        seuratObj$cluster_name <- droplevels(seuratObj$cluster_name)
        seuratObj$cluster_name_age <- paste0(seuratObj$cluster_name, "_", seuratObj$age)
        seuratObj$cluster_name_age <- factor(seuratObj$cluster_name_age)


        # create cell data object
        cds <- SeuratWrappers::as.cell_data_set(seuratObj)

        if (!useMonocleCluster) {
            fData(cds)$gene_short_name <- rownames(fData(cds))
            recreate.partition <- c(rep(1, length(cds@colData@rownames)))
            names(recreate.partition) <- cds@colData@rownames
            recreate.partition <- as.factor(recreate.partition)
            cds@clusters$UMAP$partitions <- recreate.partition

            # assign  info
            list_cluster <- seuratObj@active.ident
            cds@clusters$UMAP$clusters <- list_cluster
            cds@int_colData@listData$reducedDims$UMAP <- seuratObj@reductions$umap@cell.embeddings
            cds <- estimate_size_factors(cds)
            cds <- preprocess_cds(cds)
        } else {
            fData(cds)$gene_short_name <- rownames(fData(cds))
            cds <- estimate_size_factors(cds)
            cds <- preprocess_cds(cds, method = "PCA", num_dim = 50)
            cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")
            cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
        }

        # plot
        clusterByClusterName <- plot_cells(cds,
            color_cells_by = "cluster_name",
            label_groups_by_cluster = FALSE,
            group_label_size = 5, cell_size = 0.6
        ) + theme(legend.position = "right")
        clusterByAge <- plot_cells(cds,
            color_cells_by = "age",
            label_groups_by_cluster = FALSE,
            group_label_size = 5, cell_size = 0.6
        ) + theme(legend.position = "right")
        ggsave("before_trajectory.pdf", clusterByClusterName + clusterByAge,
            device = "pdf",
            path = figDir, width = 16, height = 8
        )

        # learn trajectory
        cds <- learn_graph(cds, use_partition = FALSE)
        clusterByClusterName <- plot_cells(cds,
            color_cells_by = "cluster_name", label_groups_by_cluster = FALSE,
            group_label_size = 4, cell_size = 0.6
        ) + theme(legend.position = "right")
        clusterByAge <- plot_cells(cds,
            color_cells_by = "age", label_groups_by_cluster = FALSE,
            group_label_size = 4, cell_size = 0.6
        ) + theme(legend.position = "right")
        ggsave("trajectory.pdf", clusterByClusterName + clusterByAge,
            device = "pdf",
            path = figDir, width = 16, height = 8
        )

        clusterByClusterNameHO <- plot_cells(subset(cds, , experiment == "Oxygen"),
            color_cells_by = "cluster_name", label_groups_by_cluster = FALSE,
            group_label_size = 4, cell_size = 0.6
        ) + theme(legend.position = "right") + ggtitle("Hyperoxia")
        clusterByClusterNameRA <- plot_cells(subset(cds, , experiment == "Room_Air"),
            color_cells_by = "cluster_name", label_groups_by_cluster = FALSE,
            group_label_size = 4, cell_size = 0.6
        ) + theme(legend.position = "right") + ggtitle("Room Air")
        ggsave("trajectoryByExp.pdf", clusterByClusterNameRA + clusterByClusterNameHO,
            device = "pdf",
            path = figDir, width = 16, height = 8
        )

        # order cell by pseudo time
        cds <- order_cells(cds, reduction_method = "UMAP", root_cells = row.names(subset(
            pData(cds),
            age == "pnd1"
        )))
        pseudo <- plot_cells(cds,
            color_cells_by = "pseudotime", label_groups_by_cluster = FALSE,
            group_label_size = 4, label_branch_points = FALSE, label_roots = FALSE,
            label_leaves = FALSE, cell_size = 0.6
        ) + theme(legend.position = "right")
        ggsave("trajectory_pseudotime.pdf", pseudo,
            device = "pdf",
            path = figDir, width = 8, height = 8
        )

        # cells ordered by monocle3 pseudotime
        cds$pseudotime <- pseudotime(cds)
        data.pseudo <- as.data.frame(colData(cds))

        pseudoByClusterAge <- ggplot(data.pseudo, aes(pseudotime, reorder(cluster_name_age, pseudotime, median),
            cluster_name_age,
            fill = cluster_name_age
        )) +
            geom_boxplot()
        pseudoByAge <- ggplot(data.pseudo, aes(pseudotime, reorder(age, pseudotime, median),
            age,
            fill = age
        )) +
            geom_boxplot()
        pseudoByCluster <- ggplot(data.pseudo, aes(pseudotime, reorder(cluster_name, pseudotime, median),
            cluster_name,
            fill = cluster_name
        )) +
            geom_boxplot()
        ggsave("boxplotByCluster_pseudotime.pdf", pseudoByClusterAge,
            device = "pdf",
            path = figDir, width = 8, height = 8
        )
        ggsave("boxplot_pseudotime.pdf", pseudoByAge + pseudoByCluster,
            device = "pdf",
            path = figDir, width = 16, height = 8
        )
        save(cds, file = file.path(saveDir, "cds.RData"))
    } else {
        cds <- loadRData(file.path(saveDir, "cds.RData"))
    }

    comparison <- "byClusterAndAge"
    dir.create(file.path(saveDir, comparison), showWarnings = FALSE)
    selected_clusters <- levels(factor(pData(cds)$cluster_name_age))
    cds <- subset(cds, , cluster_name_age %in% selected_clusters)
    PrintLog("Getting graph test result")
    if (file.exists(file.path(dataDir, "graph_test_result.csv"))) {
        pr_test_res <- read.csv(file.path(dataDir, "graph_test_result.csv"), row.names = 1)
    } else {
        pr_test_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = future::availableCores()) |>
            arrange(q_value) |>
            filter(q_value < 0.05)
        write.csv(pr_test_res, file.path(dataDir, "graph_test_result.csv"))
    }

    PrintLog("Getting pseudotime heatmap")
    pdf(paste0(figDir, "/", "heatmapPseudotime.pdf"), width = 8, height = 8)
    figs <- GetPseudotimeHeatmap(cds, pr_test_res)
    dev.off()
    ggsave("km_genes.png", path = figDir, grid.arrange(grobs = figs[-1], ncol = 2), width = 12, height = 12, units = "in")

    PrintLog("Finding gene modules")
    if (file.exists(file.path(dataDir, "gene_module.csv"))) {
        gene_module_df <- read.csv(file.path(dataDir, "gene_module.csv"), row.names = 1)
    } else {
        gene_module_df <- find_gene_modules(cds, resolution = resolution, cores = future::availableCores())
        write.csv(gene_module_df, file.path(dataDir, "gene_module.csv"))
    }

    pr_test_res %<>% rownames_to_column(var = "id")
    pr_test_res %<>% left_join(gene_module_df, by = "id")
    write.csv(pr_test_res, file = file.path(saveDir, comparison, paste0(comparison, "_genes.csv")), row.names = TRUE)

    cell_group_df <- tibble::tibble(
        cell = row.names(colData(cds)),
        cell_group = colData(cds)$cluster_name_age
    )
    agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
    module_dendro <- hclust(dist(agg_mat))
    gene_module_df$module <- factor(gene_module_df$module,
        levels = row.names(agg_mat)[module_dendro$order]
    )

    module_plots <- plot_cells(cds,
        genes = gene_module_df,
        label_cell_groups = FALSE,
        show_trajectory_graph = FALSE, cell_size = 0.6
    )
    ggsave("modules_plot.pdf", module_plots,
        device = "pdf",
        path = file.path(saveDir, comparison), width = 12, height = 10
    )

    row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
    pheatmap::pheatmap(agg_mat,
        scale = "column", clustering_method = "ward.D2",
        fontsize_row = 12, filename = file.path(saveDir, comparison, "modules_heatmap.pdf")
    )
}

#' Get a pseudotime heatmap for a CellDataSet object
#'
#' This function generates a pseudotime heatmap for a CellDataSet object.
#' It first performs a graph test on the object to identify modulated genes,
#' and then extracts the expression values for these genes and orders them by pseudotime.
#' The resulting matrix is then preprocessed and plotted as a heatmap.
#'
#' @param cds A CellDataSet object
#' @return A ComplexHeatmap plot object
#' @export
GetPseudotimeHeatmap <- function(cds, graph_test_res, morans = 0.25) {
    meta <- pData(cds)
    genes <- row.names(subset(graph_test_res, q_value == 0 & morans_I > morans))
    while (length(genes) > 200) {
        morans <- morans + 0.05
        genes <- row.names(subset(graph_test_res, q_value == 0 & morans_I > morans))
    }
    PrintLog(paste0("Number of genes selected for pseudotime heatmap: ", length(genes)))

    pt.matrix <- exprs(cds)[match(genes, rownames(rowData(cds))), order(pseudotime(cds))]
    pt.matrix <- preprocessPseudotimeHeatmap(pt.matrix)

    plots <- plotPseudotimeHeatmap(pt.matrix, "Pseudotime Heatmap", meta)
    return(plots)
}

preprocessPseudotimeHeatmap <- function(pt.matrix) {
    rowNames <- rownames(pt.matrix)
    colNames <- colnames(pt.matrix)
    pt.matrix <- t(apply(pt.matrix, 1, function(x) {
        smooth.spline(x, df = 3)$y
    }))
    pt.matrix <- t(apply(pt.matrix, 1, function(x) {
        (x - mean(x)) / sd(x)
    }))
    rownames(pt.matrix) <- rowNames
    colnames(pt.matrix) <- colNames
    pt.matrix <- pt.matrix[!rowSums(!is.finite(pt.matrix)), ]
    return(pt.matrix)
}

plotPseudotimeHeatmap <- function(pt.matrix, title, metadata, km = 4, raster = TRUE) {
    metadata <- metadata[colnames(pt.matrix), ]
    clusters <- levels(metadata$cluster_name)
    clusters_colors <- brewer.pal(9, "Set1")[1:nlevels(metadata$cluster_name)]
    colors <- list(age = c(pnd1 = "#92c44d", pnd7 = "#e67a25", pnd21 = "#b62f48"))
    colors[["cluster_name"]] <- setNames(clusters_colors, clusters)
    topAnno <- HeatmapAnnotation(
        age = metadata$age,
        cluster_name = metadata$cluster_name,
        col = colors
    )
    ht <- Heatmap(
        pt.matrix,
        name = "z-score",
        col = colorRamp2(seq(from = -2, to = 2, length = 11), rev(brewer.pal(11, "Spectral"))),
        show_row_names = TRUE,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 6),
        km = km,
        row_title_rot = 0,
        cluster_rows = ifelse(km > 1, TRUE, FALSE),
        cluster_row_slices = FALSE,
        cluster_columns = FALSE,
        column_title = title,
        top_annotation = topAnno,
        use_raster = raster
    )

    ht <- draw(ht)
    results <- list()
    plots <- list()
    plots[[1]] <- ht
    clusters <- list()

    if (km > 1) {
        pathways <- list()
        kmClusters <- row_order(ht)
        for (i in seq_along(kmClusters)) {
            plots[[length(plots) + 1]] <- GetGenesPseudotimeZscore(metadata, pt.matrix[kmClusters[[i]], ], paste0(title, " - Cluster ", i))
            pathways[[i]] <- GetEnrichRSetEnrichment(rownames(pt.matrix[kmClusters[[i]], ]), id = i)
            clusters[[i]] <- rownames(pt.matrix[kmClusters[[i]], ])
        }
        results$pathways <- do.call(rbind, pathways)
    }

    results$clusters <- clusters
    results$plots <- plots

    return(results)
}

#' Generate pseudotime heatmaps for two CellDataSet objects and save them as PDF files.
#'
#' @param cdsRA A CellDataSet object containing the RNA expression data for the reference atlas.
#' @param cdsHO A CellDataSet object containing the RNA expression data for the homology atlas.
#' @param genesRA A data frame containing the gene IDs for the reference atlas.
#' @param genesHO A data frame containing the gene IDs for the homology atlas.
#' @param cluster A string specifying the cluster name.
#' @param saveDir A string specifying the directory to save the PDF files.
#' @param append A string to append to the PDF file names.
#'
#' @return None
#'
#' @details This function generates two pseudotime heatmaps for the reference atlas and the homology atlas, respectively.
#' The heatmaps are generated based on the RNA expression data of the specified genes. The genes that are unique to each atlas are plotted separately.
#' The heatmaps are saved as PDF files in the specified directory.
#'
#' @examples
#' # Generate pseudotime heatmaps for two CellDataSet objects
#' GetPseudotimeHeatmapsByExp(cdsRA, cdsHO, genesRA, genesHO, "cluster1", "/path/to/save/directory")
#'
#' @export
GetPseudotimeHeatmapsByExp <- function(cdsRA, cdsHO, genesRA, genesHO, cluster, saveDir, append = "_1_7_21", resolutions = c(4), trajectoryPlots = FALSE) {
    metaRA <- pData(cdsRA)
    metaHO <- pData(cdsHO)

    raOnly <- setdiff(genesRA$id, genesHO$id)
    hoOnly <- setdiff(genesHO$id, genesRA$id)

    pt.matrixRA <- exprs(cdsRA)[na.omit(match(raOnly, rownames(rowData(cdsRA)))), order(pseudotime(cdsRA))]
    pt.matrixRA <- preprocessPseudotimeHeatmap(pt.matrixRA)

    pt.matrixHO <- exprs(cdsHO)[na.omit(match(hoOnly, rownames(rowData(cdsHO)))), order(pseudotime(cdsHO))]
    pt.matrixHO <- preprocessPseudotimeHeatmap(pt.matrixHO)

    for (i in resolutions) {
        pdf(paste0(saveDir, cluster, append, "/heatmapUnique_", i, ".pdf"), width = 8, height = 8)
        raGenesInRA <- plotPseudotimeHeatmap(pt.matrixRA, paste0(cluster, " - RA genes in RA"), metaRA, km = i)
        hoGenesInHO <- plotPseudotimeHeatmap(pt.matrixHO, paste0(cluster, " - HO genes in HO"), metaHO, km = i)
        dev.off()

        ggsave(paste0("RA_genes_in_RA_", i, ".pdf"),
            path = paste0(saveDir, cluster, append),
            grid.arrange(grobs = raGenesInRA$plots[-1], ncol = 2), width = 12, height = 12, units = "in"
        )
        ggsave(paste0("HO_genes_in_HO_", i, ".pdf"),
            path = paste0(saveDir, cluster, append),
            grid.arrange(grobs = hoGenesInHO$plots[-1], ncol = 2), width = 12, height = 12, units = "in"
        )

        clusterRADf <- as.data.frame(do.call(rbind, raGenesInRA$clusters))
        colnames(clusterRADf) <- names(raGenesInRA$clusters)
        write.csv(clusterRADf, paste0(saveDir, cluster, append, "/RA_genes_in_RA_", i, "_clusters.csv"))

        clusterHODf <- as.data.frame(do.call(rbind, hoGenesInHO$clusters))
        colnames(clusterHODf) <- names(hoGenesInHO$clusters)
        write.csv(clusterHODf, paste0(saveDir, cluster, append, "/HO_genes_in_HO_", i, "_clusters.csv"))

        if (!is.null(raGenesInRA$pathways)) {
            write.csv(raGenesInRA$pathways, paste0(saveDir, cluster, append, "/RA_genes_in_RA_", i, "_pathways.csv"))
        }
        if (!is.null(hoGenesInHO$pathways)) {
            write.csv(hoGenesInHO$pathways, paste0(saveDir, cluster, append, "/HO_genes_in_HO_", i, "_pathways.csv"))
        }
    }

    if (trajectoryPlots) {
        plotsRAvsHO <- GetGenePseudotimePlotByExp(cdsRA, cdsHO, raOnly)
        pdf(paste0(saveDir, cluster, append, "/RAgenes_RAvsHO.pdf"), width = 18, height = 18, onefile = TRUE)
        for (i in seq_along(plotsRAvsHO)) {
            do.call("grid.arrange", plotsRAvsHO[[i]])
        }
        dev.off()

        plotsRAvsHO <- GetGenePseudotimePlotByExp(cdsRA, cdsHO, hoOnly)
        pdf(paste0(saveDir, cluster, append, "/HOgenes_RAvsHO.pdf"), width = 18, height = 18, onefile = TRUE)
        for (i in seq_along(plotsRAvsHO)) {
            do.call("grid.arrange", plotsRAvsHO[[i]])
        }
        dev.off()
    }

    orderRA <- unlist(row_order(raGenesInRA$plots[[1]])) |> as.numeric()
    raOnly <- raOnly |>
        as.data.frame() |>
        tibble::rownames_to_column("id") |>
        dplyr::arrange(match(id, orderRA)) |>
        pull(raOnly)

    orderHO <- unlist(row_order(hoGenesInHO$plots[[1]])) |> as.numeric()
    hoOnly <- hoOnly |>
        as.data.frame() |>
        tibble::rownames_to_column("id") |>
        dplyr::arrange(match(id, orderHO)) |>
        pull(hoOnly)

    pt.matrixHORA <- exprs(cdsHO)[na.omit(match(raOnly, rownames(rowData(cdsHO)))), order(pseudotime(cdsHO))]
    pt.matrixHORA <- preprocessPseudotimeHeatmap(pt.matrixHORA)

    pt.matrixRAHO <- exprs(cdsRA)[na.omit(match(hoOnly, rownames(rowData(cdsRA)))), order(pseudotime(cdsRA))]
    pt.matrixRAHO <- preprocessPseudotimeHeatmap(pt.matrixRAHO)

    pdf(paste0(saveDir, cluster, append, "/heatmapSwapped.pdf"), width = 8, height = 8)
    hoGenesInRA <- plotPseudotimeHeatmap(pt.matrixRAHO, paste0(cluster, " - HO genes in RA"), metaRA, km = 1)
    raGenesInHO <- plotPseudotimeHeatmap(pt.matrixHORA, paste0(cluster, " - RA genes in HO"), metaHO, km = 1)
    dev.off()
}

#' Get genes pseudotime
#'
#' This function takes a monocle3 SingleCellExperiment object and a vector of gene
#' names and returns a ggplot object that shows the expression of the specified genes over pseudotime.
#'
#' @param cds A monocle3 SingleCellExperiment object
#' @param features A vector of gene names
#'
#' @return A ggplot object
#'
#' @export
GetGenesPseudotime <- function(cds, features, title) {
    cds_subset <- cds[rowData(cds)$gene_short_name %in% features, ]
    cds_subset <- cds_subset[, is.finite(colData(cds_subset)$monocle3_pseudotime)]
    cds_exprs <- SingleCellExperiment::counts(cds_subset)
    cds_exprs <- Matrix::t(Matrix::t(cds_exprs) / size_factors(cds_subset))
    cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
    colnames(cds_exprs) <- c("f_id", "Cell", "expression")
    cds_colData <- colData(cds_subset) |>
        as.data.frame() |>
        tibble::rownames_to_column("Cell")
    cds_rowData <- rowData(cds_subset) |>
        as.data.frame() |>
        tibble::rownames_to_column("f_id")
    cds_exprs <- merge(cds_exprs, cds_rowData, by = "f_id")
    cds_exprs <- merge(cds_exprs, cds_colData, by = "Cell")
    cds_exprs$adjusted_expression <- cds_exprs$expression
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
        cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
        cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    } else {
        cds_exprs$feature_label <- cds_exprs$f_id
    }
    cds_exprs$f_id <- as.character(cds_exprs$f_id)
    cds_exprs$feature_label <- factor(cds_exprs$feature_label)
    new_data <- as.data.frame(colData(cds_subset))
    new_data$Size_Factor <- 1
    model_tbl <- fit_models(cds_subset, model_formula_str = "~ splines::ns(monocle3_pseudotime, df=3)")
    model_expectation <- model_predictions(model_tbl, new_data = new_data)
    colnames(model_expectation) <- colnames(cds_subset)
    expectation <- plyr::ddply(
        cds_exprs, plyr::.(f_id, Cell),
        function(x) {
            data.frame(expectation = model_expectation[
                x$f_id,
                x$Cell
            ])
        }
    )
    cds_exprs <- merge(cds_exprs, expectation)

    q <- ggplot(aes(monocle3_pseudotime, expression), data = cds_exprs) +
        geom_line(aes(x = monocle3_pseudotime, y = expectation, group = feature_label), data = cds_exprs) +
        ggtitle(title)
    return(q)
}

#' Get gene expression plot over pseudotime
#'
#' This function generates a plot of gene expression over pseudotime for a given cell dataset and data matrix.
#'
#' @param metadata A monocle3 colData object containing the pseudotime information as monocle3_pseudotime.
#' @param dataMatrix A data matrix containing the gene expression data.
#'
#' @return A ggplot object containing the gene expression plot.
#'
#' @export
GetGenesPseudotimeZscore <- function(metadata, dataMatrix, title) {
    pseudotime <- metadata[colnames(dataMatrix), "monocle3_pseudotime", drop = FALSE] |>
        as.data.frame() |>
        dplyr::rename(pseudotime = monocle3_pseudotime)
    dataMatrix <- dataMatrix %>%
        t() %>%
        as.data.frame() %>%
        merge(pseudotime, by = "row.names") %>%
        tibble::column_to_rownames("Row.names") %>%
        dplyr::mutate(mean = rowMeans(dplyr::select(., -pseudotime), na.rm = TRUE)) %>%
        pivot_longer(cols = -pseudotime, names_to = "gene_name", values_to = "expression") %>%
        as.data.frame()

    # "loess"
    dataMatrix$type <- ifelse(dataMatrix$gene_name == "mean", "mean", "gene")
    dataMatrix$size <- ifelse(dataMatrix$type == "gene", 0.5, 2)
    q <- ggplot(data = dataMatrix, aes(x = pseudotime, y = expression, group = gene_name, color = type)) +
        geom_smooth(se = FALSE, aes(size = size), method = "loess") +
        theme_bw() +
        theme(
            panel.border = element_rect(fill = NA, color = "black"),
            panel.grid = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        ) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_size_identity() +
        scale_color_manual(values = c("gray", "red")) +
        ggtitle(title)

    return(q)
}

#' Get gene pseudotime plots for two SingleCellExperiment objects
#'
#' This function takes in two SingleCellExperiment objects (cdsRA and cdsHO)
#' and a vector of gene names (geneList), and returns a list of plots showing
#' the expression of the genes in pseudotime for each dataset.
#'
#' @param cdsRA A SingleCellExperiment object containing the RA dataset.
#' @param cdsHO A SingleCellExperiment object containing the HO dataset.
#' @param geneList A vector of gene names.
#'
#' @return A list of plots showing the expression of the genes in pseudotime for each dataset.
#'
#' @examples
#' cdsRA <- SingleCellExperiment(assays = list(counts = matrix(rnbinom(1000, 100, 0.1), ncol = 10)))
#' cdsHO <- SingleCellExperiment(assays = list(counts = matrix(rnbinom(1000, 100, 0.1), ncol = 10)))
#' geneList <- c("gene1", "gene2", "gene3", "gene4", "gene5", "gene6")
#' GetGenePseudotimePlotByExp(cdsRA, cdsHO, geneList)
#'
#' @seealso \code{\link{plot_genes_in_pseudotime}}
GetGenePseudotimePlotByExp <- function(cdsRA, cdsHO, geneList) {
    plotList <- list()
    for (i in seq(1, length(geneList), 6)) {
        genes <- geneList[i:min(i + 6 - 1, length(geneList))]
        raPlot <- plot_genes_in_pseudotime(cdsRA[rowData(cdsRA)$gene_short_name %in% genes, ],
            color_cells_by = "age",
            min_expr = 0.5
        )
        hoPlot <- plot_genes_in_pseudotime(cdsHO[rowData(cdsHO)$gene_short_name %in% genes, ],
            color_cells_by = "age",
            min_expr = 0.5
        )
        plotList[[length(plotList) + 1]] <- grid.arrange(raPlot, hoPlot, ncol = 2)
    }

    finalList <- list()
    for (j in seq(1, length(plotList), 3)) {
        finalList[[length(finalList) + 1]] <- grid.arrange(grobs = plotList[j:min(j + 3 - 1, length(plotList))], ncol = 3, nrow = 1)
    }

    return(finalList)
}
