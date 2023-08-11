#' Get clusters for a Seurat sample using PCA, UMAP, and t-SNE
#'
#' @param seuratSample A Seurat object containing gene expression data for a specific sample
#' @param TSNE A logical value indicating whether to run t-SNE
#' @param verbose A logical value indicating whether to print verbose output
#' @param clustersFrom A numeric value specifying the minimum resolution for clustering
#' @param clustersTo A numeric value specifying the maximum resolution for clustering
#' @param clustersBy A numeric value specifying the resolution increment for clustering
#' @param stepName A character string specifying the name of the analysis step
#' @param age A character string specifying the age of the sample
#'
#' @return A list containing the Seurat object with clusters and a ggplot object
#'
#' @export
GetSeuratClusters <- function(seuratSample, TSNE = TRUE, verbose = FALSE,
                              clustersFrom = 0.4, clustersTo = 1.0, clustersBy = 0.1,
                              stepName = "") {
    seuratSample <- Seurat::RunPCA(seuratSample, verbose = verbose)
    dimsCalc <- GetSignificantPC(seuratSample, stepName, verbose = verbose)
    seuratSample <- Seurat::RunUMAP(seuratSample, dims = 1:dimsCalc, verbose = verbose)
    if (TSNE) {
        seuratSample <- Seurat::RunTSNE(seuratSample, dims.use = 1:dimsCalc, verbose = verbose)
    }
    seuratSample <- Seurat::FindNeighbors(seuratSample, dims = 1:dimsCalc, reduction = "pca", verbose = verbose)
    seuratSample <- Seurat::FindClusters(seuratSample, resolution = seq(from = clustersFrom, to = clustersTo, by = clustersBy), verbose = verbose)
    return(seuratSample)
}

#' Get the number of significant principal components for a Seurat sample
#'
#' @param seuratSample A Seurat object containing gene expression data for a specific sample
#' @param stepName A character string specifying the name of the analysis step
#' @param age A character string specifying the age of the sample
#' @param verbose A logical value indicating whether to print verbose output
#'
#' @return A list containing the number of significant principal components and a ggplot object
#'
#' @export
GetSignificantPC <- function(seuratSample, stepName = "", verbose = FALSE) {
    stdv <- seuratSample[["pca"]]@stdev
    percent.stdv <- (stdv / sum(stdv)) * 100
    cumulative <- cumsum(percent.stdv)
    co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
    co2 <- sort(which((percent.stdv[seq_along(percent.stdv) - 1] - percent.stdv[2:length(percent.stdv)]) > 0.1), decreasing = TRUE)[1] + 1
    min.pc <- min(co1, co2)
    p <- Seurat::ElbowPlot(seuratSample, ndims = 40) + ggplot2::geom_vline(xintercept = min.pc)
    if (stepName != "" && age != "") {
        ggplot2::ggsave(
            filename = paste0("elbow_", stepName, ".pdf"),
            device = "pdf",
            plot = p,
            path = file.path("./figures", age, "QC", seuratSample@project.name),
            width = 8, height = 8, units = "in"
        )
    }
    return(min.pc)
}

#' Get a subset of cells from a Seurat object based on cell type
#'
#' This function takes a Seurat object and a cell type name, and returns a subset of the Seurat object containing only
#' cells of the specified cell type. The function performs SCTransform normalization on the subset and runs the Seurat
#' clustering algorithm on the normalized data. The resulting Seurat object has cluster identities assigned to it and is returned.
#'
#' @param seuratObj A Seurat object containing the cells to subset.
#' @param cellTypeName A character string specifying the name of the cell type to subset.
#' @param verbose A logical value indicating whether to print verbose output.
#'
#' @return A Seurat object containing only cells of the specified cell type.
#'
#' @export
GetCellTypeSubset <- function(seuratObj, cellTypeName, verbose = FALSE) {
    seuratSubset <- subset(seuratObj, subset = celltype == cellTypeName)
    seuratSubset <- SCTransform(seuratSubset, vst.flavor = "v2", return.only.var.genes = FALSE, verbose = verbose)
    seuratSubset <- GetSeuratClusters(seuratSubset, stepName = paste0(cellTypeName, "_subset"))

    Idents(seuratSubset) <- "cluster_name"
    return(seuratSubset)
}

#' Get the number of cells in each cluster for different factors
#'
#' @param seuratObj A Seurat object containing gene expression data
#' @param extraFactor An optional character string specifying an additional factor to group by
#'
#' @return A data frame containing the number of cells in each cluster for each factor
#'
#' @export
GetCellNumber <- function(seuratObj, extraFactor = NULL) {
    seuratObj@meta.data$cluster_names <- seuratObj@active.ident
    md <- seuratObj@meta.data |> data.table::as.data.table()
    cell_num <- md[, .N, by = c("sample", "cluster_names")] %>%
        data.table::dcast(., sample ~ cluster_names, value.var = "N") %>%
        t() %>%
        base::as.data.frame()
    colnames(cell_num) <- cell_num[1, ]
    cell_num <- cell_num[-1, , drop = FALSE]

    cell_num3 <- md[, .N, by = c("sex", "cluster_names")] %>%
        data.table::dcast(., sex ~ cluster_names, value.var = "N") %>%
        t() %>%
        base::as.data.frame()
    colnames(cell_num3) <- cell_num3[1, ]
    cell_num3 <- cell_num3[-1, , drop = FALSE]
    cell_num3 <- cell_num3 %>% tibble::rownames_to_column(var = "cluster")

    cell_num4 <- md[, .N, by = c("experiment", "cluster_names")] %>%
        data.table::dcast(., experiment ~ cluster_names, value.var = "N") %>%
        t() %>%
        base::as.data.frame()
    colnames(cell_num4) <- cell_num4[1, ]
    cell_num4 <- cell_num4[-1, , drop = FALSE]
    cell_num4 <- cell_num4 %>% tibble::rownames_to_column(var = "cluster")

    cell_num5 <- md[, .N, by = c("age", "cluster_names")] %>%
        data.table::dcast(., age ~ cluster_names, value.var = "N") %>%
        t() %>%
        base::as.data.frame()
    colnames(cell_num5) <- cell_num5[1, ]
    cell_num5 <- cell_num5[-1, , drop = FALSE]
    cell_num5 <- cell_num5 %>% tibble::rownames_to_column(var = "cluster")

    a <- as.array(table(Seurat::Idents(seuratObj)))
    cell_num$integrated <- a
    cell_num <- cell_num |>
        tibble::rownames_to_column(var = "cluster") |>
        dplyr::left_join(cell_num3, by = "cluster") |>
        dplyr::left_join(cell_num4, by = "cluster") |>
        dplyr::left_join(cell_num5, by = "cluster")

    extraLevels <- c()
    if (!is.null(extraFactor)) {
        if (extraFactor %in% colnames(seuratObj@meta.data)) {
            cell_num6 <- md[, .N, by = c(extraFactor, "cluster_names")] %>%
                data.table::dcast(., base::get(extraFactor) ~ cluster_names, value.var = "N") %>%
                t() %>%
                base::as.data.frame()
            cell_num6[1, ] <- ifelse(is.na(cell_num6[1, ]), paste0(extraFactor, "_None"), paste0(extraFactor, "_", cell_num6[1, ]))
            colnames(cell_num6) <- cell_num6[1, ]
            cell_num6 <- cell_num6[-1, , drop = FALSE]
            cell_num6 <- cell_num6 %>% tibble::rownames_to_column(var = "cluster")
            cell_num <- cell_num |> dplyr::left_join(cell_num6, by = "cluster")
            extraLevels <- levels(factor(seuratObj@meta.data[[extraFactor]]))
            extraLevels <- paste0(extraFactor, "_", extraLevels)
        } else {
            print("Error: Specified extraFactor is not present in metadata")
        }
    }

    cell_num <- cell_num %>% tibble::column_to_rownames(var = "cluster")
    cell_num <- cell_num %>% dplyr::select(c(
        "integrated", levels(factor(seuratObj$experiment)),
        levels(factor(seuratObj$sex)), levels(factor(seuratObj$age)),
        levels(factor(seuratObj$sample)), extraLevels
    ))
    cell_num[] <- sapply(cell_num, as.numeric)
    cell_num <- cell_num %>% replace(is.na(.), 0)
    return(cell_num)
}

#' Get x-axis tick values from a ggplot object
#'
#' This function takes a ggplot object as input and returns the x-axis tick values as a numeric vector.
#'
#' @param plot a ggplot object
#' @return A numeric vector containing the x-axis tick values of the plot.
#' @examples
#' plot <- ggplot(mtcars, aes(x = mpg, y = wt)) +
#'     geom_point()
#' get_x_ticks(plot)
#' @export
get_x_ticks <- function(plot) {
    plot <- ggplot2::ggplot_build(plot)
    return(na.omit(plot$layout$panel_params[[1]]$x$breaks))
}

#' Get y-axis tick values from a ggplot object
#'
#' This function takes a ggplot object as input and returns the y-axis tick values as a numeric vector.
#'
#' @param plot a ggplot object
#' @return A numeric vector containing the y-axis tick values of the plot.
#' @examples
#' plot <- ggplot(mtcars, aes(x = mpg, y = wt)) +
#'     geom_point()
#' get_y_ticks(plot)
#' @export
get_y_ticks <- function(plot) {
    plot <- ggplot2::ggplot_build(plot)
    return(na.omit(plot$layout$panel_params[[1]]$y$breaks))
}

#' Generate a vector of colours suitable for use in ggplot2 plots
#'
#' This function generates a vector of n colours suitable for use in ggplot2 plots.
#'
#' @param n the number of colours to generate (default is 6)
#' @param h a numeric vector of length 2 specifying the range of hues to use (default is c(0, 360) + 15)
#' @return A vector of n colours in hexadecimal format.
#' @examples
#' ggplotColours(3, c(0, 240))
#' @export
ggplotColours <- function(n = 6, h = c(0, 360) + 15) {
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360 / n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

#' Shorten given cluster name
#'
#' @param clusterName string to shorten
#' @return String
#' @examples
#' shortenClusterName("Oligodendrocyte precursor cell") # returns "Oli_pre_cel"
#' @export
shortenClusterName <- function(clusterName) {
    return(paste(clusterName |>
        strsplit(split = " ") |>
        unlist() |>
        purrr::map_if(~ !grepl("\\+", .x),
            \(x) substr(x, start = 1, stop = 3),
            .else = as.character()
        ) |>
        unlist(), collapse = " "))
}

#' Remove special characters and spaces form cluster names
#'
#' @param clusterName string to simplify
#' @return String
#' @examples
#' simplifyClusterName("Oligodendrocyte precursor cell") # returns "Oligodendrocyte_precursor_cell"
#' @export
SimplifyClusterName <- function(clusterName) {
    clusterNameUnderscore <- gsub(" ", "_", clusterName, fixed = TRUE)
    return(gsub("+", "Pos", clusterNameUnderscore, fixed = TRUE))
}
