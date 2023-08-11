#' Run Cellchat on multiple datasets
#'
#' This function generates cellchat objects and figures
#'
#' @param seuratObj Seurat object with hyperoxia and Room air cells
#' @param cellType Celltype of Seurat cluster being analyzed
#' @param age Age of sample being analyzed
#' @param byExp Wheter we need to create two cellchat objects RA & HO
#' @return Null
#' @export
RunCellChat <- function(seuratObj, cellType, age, byExp = TRUE, projectDir) {
    # Set options to run cellchat
    options(stringsAsFactors = FALSE)
    options(parallelly.makeNodePSOCK.validate = (packageVersion("parallelly") > "1.31.1"))

    # prepare cluster names for analysis
    seuratObj$cluster_names <- seuratObj@active.ident

    # start excel workbook
    wb <- openxlsx::createWorkbook()
    headerStyle <- openxlsx::createStyle(
        textDecoration = "Bold",
        fontSize = 14, fontColour = "#FFFFFF", halign = "center",
        fgFill = "#4F81BD", border = "TopBottom", borderColour = "#4F81BD"
    )

    # create directories to store results
    dir.create(file.path(projectDir, "data", age, "byCellType", cellType), recursive = TRUE)
    dir.create(file.path(projectDir, "figures", age, "byCellType", cellType), recursive = TRUE)
    if (byExp) {
        # load cellchat objects
        cellchatOxy <- GetCellChatObject(seuratObj, cellType, age, "cellchatOxy.RData", "Oxygen", projectDir = projectDir)
        cellchatRoom <- GetCellChatObject(seuratObj, cellType, age, "cellchatRoom.RData", "Room_Air", projectDir = projectDir)

        savePath <- file.path(projectDir, "figures", age, "byCellType", cellType, "cellchat/byExp/")
        dir.create(savePath, recursive = TRUE)
        dir.create(paste0(savePath, "netVisual_chord/"))

        # merge objects
        object.list <- list(roomair = cellchatRoom, oxygen = cellchatOxy)
        cellchat <- CellChat::mergeCellChat(object.list, add.names = names(object.list))
        cellchat <- CellChat::liftCellChat(cellchat, union(levels(cellchatRoom@idents), levels(cellchatOxy@idents)))

        # save communication data as excel spreadsheet
        room <- subsetCommunication(cellchat)$roomair |> as.data.frame()
        oxy <- subsetCommunication(cellchat)$oxygen |> as.data.frame()
        openxlsx::addWorksheet(wb, sheetName = paste(cellType, "Room Air"), gridLines = TRUE, zoom = 180)
        openxlsx::writeData(wb,
            sheet = paste(cellType, "Room Air"), x = room,
            headerStyle = headerStyle, rowNames = FALSE
        )
        openxlsx::addWorksheet(wb, sheetName = paste(cellType, "Oxygen"), gridLines = TRUE, zoom = 180)
        openxlsx::writeData(wb,
            sheet = paste(cellType, "Oxygen"), x = oxy,
            headerStyle = headerStyle, rowNames = FALSE
        )
        openxlsx::saveWorkbook(wb,
            file = paste0(savePath, cellType, "_inferred_communication.xlsx"),
            overwrite = TRUE
        )

        # generate plots
        gg1 <- CellChat::compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2))
        gg2 <- CellChat::compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "weight")
        ggplot2::ggsave("compareInteractions.pdf", gg1 + gg2,
            device = "pdf", height = 12, width = 12,
            path = savePath
        )

        par(mfrow = c(1, 2), xpd = TRUE)
        pdf(paste0(savePath, "newVisual_diffInteraction.pdf"), width = 8, heigh = 8)
        CellChat::netVisual_diffInteraction(cellchat, weight.scale = TRUE)
        CellChat::netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
        dev.off()

        gg1 <- CellChat::netVisual_heatmap(cellchat)
        gg2 <- CellChat::netVisual_heatmap(cellchat, measure = "weight")
        pdf(paste0(savePath, "netVisual_heatmap.pdf"), width = 8, heigh = 8)
        ComplexHeatmap::draw(gg1)
        ComplexHeatmap::draw(gg2)
        dev.off()

        weight.max <- CellChat::getMaxWeight(object.list, attribute = c("idents", "weight"))
        par(mfrow = c(1, 2), xpd = TRUE)
        for (i in seq_along(object.list)) {
            pdf(paste0(savePath, names(object.list)[i], "_netVisual_circle.pdf"), width = 8, heigh = 8)
            CellChat::netVisual_circle(object.list[[i]]@net$weight,
                weight.scale = TRUE, label.edge = FALSE, edge.weight.max = weight.max[2],
                edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i])
            )
            dev.off()
        }

        num.link <- sapply(object.list, function(x) {
            rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
        }) |> unlist()
        weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
        gg <- list()
        for (i in seq_along(object.list)) {
            gg[[i]] <- CellChat::netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
        }
        ggsave("netAnalysis_signalingRole_scatter.pdf", patchwork::wrap_plots(plots = gg),
            device = "pdf", height = 8, width = 12,
            path = savePath
        )

        cellchat <- CellChat::computeNetSimilarityPairwise(cellchat, type = "functional")
        #> Compute signaling network similarity for datasets 1 2
        cellchat <- CellChat::netEmbedding(cellchat, type = "functional")
        #> Manifold learning of the signaling networks for datasets 1 2
        cellchat <- CellChat::netClustering(cellchat, type = "functional", do.parallel = FALSE)
        #> Classification learning of the signaling networks for datasets 1 2
        # Visualization in 2D-space
        gg1 <- CellChat::netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 2, do.label = TRUE, top.label = .8)
        ggplot2::ggsave("netVisual_embeddingPairwise_functional.pdf", gg1,
            device = "pdf", height = 8, width = 12,
            path = savePath
        )

        cellchat <- CellChat::computeNetSimilarityPairwise(cellchat, type = "structural")
        #> Compute signaling network similarity for datasets 1 2
        cellchat <- CellChat::netEmbedding(cellchat, type = "structural")
        #> Manifold learning of the signaling networks for datasets 1 2
        cellchat <- CellChat::netClustering(cellchat, type = "structural", do.parallel = FALSE)
        #> Classification learning of the signaling networks for datasets 1 2
        # Visualization in 2D-space
        gg1 <- CellChat::netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 2, do.label = TRUE, top.label = .8)
        ggplot2::ggsave("netVisual_embeddingPairwise_structural.pdf", gg1,
            device = "pdf", height = 8, width = 12,
            path = savePath
        )

        gg1 <- CellChat::netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
        ggplot2::ggsave("netVisual_embeddingPairwiseZoomIn.pdf", gg1,
            device = "pdf", height = 14, width = 24,
            path = savePath
        )

        gg1 <- CellChat::rankSimilarity(cellchat, type = "functional")
        ggplot2::ggsave("rankSimilarity.pdf", gg1,
            device = "pdf", height = 8, width = 12,
            path = savePath
        )

        gg1 <- CellChat::rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)
        gg2 <- CellChat::rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE)
        ggplot2::ggsave("rankNet.pdf", gg1 + gg2,
            device = "pdf", height = 12, width = 12,
            path = savePath
        )

        i <- 1
        # combining all the identified signaling pathways from different datasets
        pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i + 1]]@netP$pathways)
        ht1 <- CellChat::netAnalysis_signalingRole_heatmap(object.list[[i]],
            pattern = "outgoing",
            signaling = pathway.union, title = names(object.list)[i], width = 8, height = 12
        )
        ht2 <- CellChat::netAnalysis_signalingRole_heatmap(object.list[[i + 1]],
            pattern = "outgoing",
            signaling = pathway.union, title = names(object.list)[i + 1], width = 8, height = 12
        )
        pdf(paste0(savePath, "netAnalysis_signalingRole_heatmap_outgoing.pdf"), width = 12, heigh = 8)
        print(ht1 + ht2, ht_gap = unit(0.5, "cm"))
        dev.off()

        ht1 <- CellChat::netAnalysis_signalingRole_heatmap(object.list[[i]],
            pattern = "incoming",
            signaling = pathway.union, title = names(object.list)[i], width = 8, height = 12, color.heatmap = "GnBu"
        )
        ht2 <- CellChat::netAnalysis_signalingRole_heatmap(object.list[[i + 1]],
            pattern = "incoming",
            signaling = pathway.union, title = names(object.list)[i + 1], width = 8, height = 12, color.heatmap = "GnBu"
        )
        pdf(paste0(savePath, "netAnalysis_signalingRole_heatmap_incoming.pdf"), width = 12, heigh = 8)
        print(ht1 + ht2, ht_gap = unit(0.5, "cm"))
        dev.off()

        ht1 <- CellChat::netAnalysis_signalingRole_heatmap(object.list[[i]],
            pattern = "all",
            signaling = pathway.union, title = names(object.list)[i], width = 8, height = 12, color.heatmap = "OrRd"
        )
        ht2 <- CellChat::netAnalysis_signalingRole_heatmap(object.list[[i + 1]],
            pattern = "all",
            signaling = pathway.union, title = names(object.list)[i + 1], width = 8, height = 12, color.heatmap = "OrRd"
        )
        pdf(paste0(savePath, "netAnalysis_signalingRole_heatmap_all.pdf"), width = 12, heigh = 8)
        print(ht1 + ht2, ht_gap = unit(0.5, "cm"))
        dev.off()

        allClusters <- levels(object.list[[1]]@idents)
        for (i in seq_along(object.list)) {
            for (clusterName in allClusters) {
                tryCatch(
                    {
                        pdf(paste0(savePath, "netVisual_chord/", names(object.list)[i], "_", clusterName, "_netVisual_chord_gene.pdf"), width = 12, heigh = 12)
                        netVisual_chord_gene(object.list[[i]],
                            sources.use = clusterName, targets.use = allClusters[allClusters != clusterName],
                            lab.cex = 0.5, title.name = paste0("Signaling from ", clusterName, " - ", names(object.list)[i]), small.gap = 1, big.gap = 1
                        )
                        dev.off()
                    },
                    error = function(e) {
                        cat("ERROR :", conditionMessage(e), "\n")
                    }
                )
            }
        }

        # Get plot to select cluster number based in similarity
        # Create a cluster with two workers
        cl <- parallel::makeCluster(2)
        # Create a list of the two patterns you want to use
        patterns <- c("outgoing", "incoming")
        # Use mclapply to run the function on both patterns in parallel
        plots <- parallel::mclapply(patterns, SelectKParallel, cellchat = cellchatRoom, mc.cores = 2)
        # Assign the results to the appropriate variables
        outPlotRoom <- plots[[1]]
        inPlotRoom <- plots[[2]]
        # Stop the cluster
        parallel::stopCluster(cl)
        outkRoom <- GetBestK(outPlotRoom)
        inkRoom <- GetBestK(inPlotRoom)
        ggplot2::ggsave("roomK.pdf", outPlotRoom + inPlotRoom, device = "pdf", path = savePath)
        cellchatRoom <- GetPatternsAndSimilarity(cellchatRoom, outkRoom, inkRoom)

        save(cellchatRoom, file = paste0("./data/", age, "/byCellType/", cellType, "/cellchat_room_computed.RData"))
        saveRDS(cellchatRoom, paste0("./data/", age, "/byCellType/", cellType, "/cellchat_room_computed.rds"))

        gg1 <- netAnalysis_river(cellchatRoom, pattern = "outgoing", font.size = 4, font.size.title = 14)
        gg2 <- netAnalysis_river(cellchatRoom, pattern = "incoming", font.size = 4, font.size.title = 14)
        ggsave("room_out_communication.pdf", gg1, device = "pdf", path = savePath)
        ggsave("room_in_communication.pdf", gg2, device = "pdf", path = savePath)

        # Get plot to select cluster number based in similarity
        # Create a cluster with two workers
        cl <- parallel::makeCluster(2)
        # Create a list of the two patterns you want to use
        patterns <- c("outgoing", "incoming")
        # Use mclapply to run the function on both patterns in parallel
        plots <- parallel::mclapply(patterns, SelectKParallel, cellchat = cellchatOxy, mc.cores = 2)
        # Assign the results to the appropriate variables
        outPlotOxy <- plots[[1]]
        inPlotOxy <- plots[[2]]
        # Stop the cluster
        parallel::stopCluster(cl)
        outkOxy <- GetBestK(outPlotOxy)
        inkOxy <- GetBestK(inPlotOxy)
        ggplot2::ggsave("oxyK.pdf", outPlotOxy + inPlotOxy, device = "pdf", path = savePath)
        cellchatOxy <- GetPatternsAndSimilarity(cellchatOxy, outkOxy, inkOxy)

        save(cellchatOxy, file = paste0("./data/", age, "/byCellType/", cellType, "/cellchat_oxy_computed.RData"))
        saveRDS(cellchatOxy, paste0("./data/", age, "/byCellType/", cellType, "/cellchat_oxy_computed.rds"))

        gg1 <- netAnalysis_river(cellchatOxy, pattern = "outgoing", font.size = 4, font.size.title = 14)
        gg2 <- netAnalysis_river(cellchatOxy, pattern = "incoming", font.size = 4, font.size.title = 14)
        ggsave("oxy_out_communication.pdf", gg1, device = "pdf", path = savePath)
        ggsave("oxy_in_communication.pdf", gg2, device = "pdf", path = savePath)
    } else {
        # load cellchat objects
        cellchat <- GetCellChatObject(seuratObj, cellType, age, projectDir = projectDir)

        savePath <- file.path(projectDir, "figures", age, "byCellType", cellType, "cellchat/")
        dir.create(savePath)
        dir.create(paste0(savePath, "netVisual_chord/"))

        # save communication data as excel spreadsheet
        openxlsx::addWorksheet(wb, sheetName = cellType, gridLines = TRUE, zoom = 180)
        openxlsx::writeData(wb,
            sheet = cellType, x = subsetCommunication(cellchat),
            headerStyle = headerStyle, rowNames = FALSE
        )
        openxlsx::saveWorkbook(wb,
            file = paste0(savePath, cellType, "_inferred_communication.xlsx"),
            overwrite = TRUE
        )

        # Get plot to select cluster number based in similarity
        # Create a cluster with two workers
        cl <- parallel::makeCluster(2)
        # Create a list of the two patterns you want to use
        patterns <- c("outgoing", "incoming")
        # Use mclapply to run the function on both patterns in parallel
        plots <- mclapply(patterns, SelectKParallel, cellchat = cellchat, mc.cores = 2)
        # Assign the results to the appropriate variables
        out_plot <- plots[[1]]
        in_plot <- plots[[2]]
        # Stop the cluster
        parallel::stopCluster(cl)
        outk <- GetBestK(out_plot)
        ink <- GetBestK(in_plot)
        ggsave("ccK.pdf", out_plot + in_plot, device = "pdf", path = savePath, width = 12, heigh = 8)

        # Generate plots
        groupSize <- as.numeric(table(cellchat@idents))
        par(mfrow = c(1, 2), xpd = TRUE)
        pdf(paste0(savePath, "netVisual_circle.pdf"), width = 8, heigh = 8)
        CellChat::netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
        CellChat::netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")
        dev.off()

        ht1 <- CellChat::netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
        ht2 <- CellChat::netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
        pdf(paste0(savePath, "netAnalysis_signalingRole_heatmap.pdf"), width = 12, heigh = 8)
        ComplexHeatmap::draw(ht1)
        ComplexHeatmap::draw(ht2)
        dev.off()

        cellchat <- GetPatternsAndSimilarity(cellchat, outk, ink)
        save(cellchat, file = paste0("./data/", age, "/byCellType/", cellType, "/cellchat_computed.RData"))
        saveRDS(cellchat, paste0("./data/", age, "/byCellType/", cellType, "/cellchat_computed.rds"))

        gg1 <- netAnalysis_river(cellchat, pattern = "outgoing", font.size = 4, font.size.title = 14)
        gg2 <- netAnalysis_river(cellchat, pattern = "incoming", font.size = 4, font.size.title = 14)
        ggsave("out_communication.pdf", gg1, device = "pdf", path = savePath)
        ggsave("in_communication.pdf", gg2, device = "pdf", path = savePath)
    }
}

#' Get or create a CellChat object
#'
#' This function checks if a CellChat object exists in the specified file path.
#' If it does, it loads the object and returns it.
#' If it doesn't, it runs the preprocessing pipeline to create a new CellChat object and saves it to the specified file path.
#'
#' @param seuratObj Seurat object. If subset is set hyperoxia and Room air cells are required.
#' @param cellType Type of cell to analyze
#' @param age Age of the sample
#' @param filename Name of the file to save the CellChat object to
#' @param subset An optional subset argument to filter the Seurat object before running the preprocessing pipeline. Uses experiment property.
#' @return CellChat object
#' @export
GetCellChatObject <- function(seuratObj, cellType, age, filename = "cellchat.RData", subset, projectDir) {
    if (!file.exists(paste0("./data/", age, "/byCellType/", cellType, "/", filename))) {
        # If the file doesn't exist, run the preprocessing pipeline to create a new cellchat object and save it to the specified file path.
        if (!missing(subset)) {
            seuratObj <- subset(seuratObj, subset = experiment == subset)
        }
        seuratObj <- RemoveSmallClusters(seuratObj)
        cellchat <- PreprocessCellChat(seuratObj)
        save(cellchat, file = file.path(projectDir, "data", age, "byCellType", cellType, filename))
    } else {
        # load cellchat objects
        load(file.path(projectDir, "data", age, "byCellType", cellType, filename))
    }
    return(cellchat)
}

#' Pre process Cellchat object
#'
#' This function runs all the necessary functions to start analysing a cellchat object
#'
#' @param seuratObj Seurat object with hyperoxia and Room air cells
#' @param experiment name of analysis used to print log messages
#' @return CellChat object
#' @export
PreprocessCellChat <- function(seuratObj, experiment = "CellChat") {
    data.input <- Seurat::GetAssayData(seuratObj, assay = "SCT", slot = "data")
    # ignore empty levels (a cluster with no cells)
    seuratObj$cluster_names <- droplevels(seuratObj$cluster_names, exclude = setdiff(levels(seuratObj$cluster_names), unique(seuratObj$cluster_names)))
    # create cellchat object
    cellchatObj <- CellChat::createCellChat(object = data.input, meta = seuratObj@meta.data, group.by = "cluster_names")
    # set cellchat database to mouse
    cellchatObj@DB <- CellChat::CellChatDB.mouse

    # Run preprocessing steps
    cellchatObj <- CellChat::subsetData(cellchatObj) # This step is necessary even if using the whole database
    cellchatObj <- CellChat::identifyOverExpressedGenes(cellchatObj)
    cellchatObj <- CellChat::identifyOverExpressedInteractions(cellchatObj)
    cellchatObj <- CellChat::projectData(cellchatObj, PPI.mouse)
    cellchatObj <- CellChat::computeCommunProb(cellchatObj)
    cellchatObj <- CellChat::filterCommunication(cellchatObj, min.cells = 100)
    cellchatObj <- CellChat::computeCommunProbPathway(cellchatObj)
    cellchatObj <- CellChat::aggregateNet(cellchatObj)
    cellchatObj <- CellChat::netAnalysis_computeCentrality(cellchatObj)
    return(cellchatObj)
}

#' Set communication patterns on CellChat
#'
#' Set cluster number for communication patterns and run necessary functions
#' for further analysis
#'
#' @param seuratObj Seurat object with hyperoxia and Room air cells
#' @param experiment name of analysis used to print log messages
#' @return CellChat object
#' @export
GetPatternsAndSimilarity <- function(cellChatObject, pattern_out, pattern_in, umap.method = "uwot") {
    # Set number of clusters
    cellchat <- CellChat::identifyCommunicationPatterns(cellChatObject, pattern = "outgoing", k = pattern_out)
    cellchat <- CellChat::identifyCommunicationPatterns(cellchat, pattern = "incoming", k = pattern_in)

    cellchat <- CellChat::computeNetSimilarity(cellchat, type = "functional")
    cellchat <- CellChat::netEmbedding(cellchat, type = "functional", umap.method = umap.method)
    cellchat <- CellChat::netClustering(cellchat, type = "functional", do.parallel = FALSE)

    cellchat <- CellChat::computeNetSimilarity(cellchat, type = "structural")
    cellchat <- CellChat::netEmbedding(cellchat, type = "structural", umap.method = umap.method)
    cellchat <- CellChat::netClustering(cellchat, type = "structural", do.parallel = FALSE)
    return(cellchat)
}

#' Get Male vs Female Comparison
#'
#' This function runs all the necessary functions to start analysing a cellchat female object
#' vs a male cellchat object
#'
#' @param cluster_name name of cluster to be compared
#' @param maleCellChat male cellchat object
#' @param femaleCellChat female cellchat object
#' @param age age of the mice
#' @return NULL
#' @export
GetCellChatComparisonBySex <- function(cluster_name, maleCellChat, femaleCellChat, age) {
    object.list <- list(female = femaleCellChat, male = maleCellChat)
    cellchat <- CellChat::mergeCellChat(object.list, add.names = names(object.list))
    cellchat <- CellChat::liftCellChat(cellchat, union(levels(femaleCellChat@idents), levels(maleCellChat@idents)))

    gg1 <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)
    gg2 <- rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE)
    ggsave("rankNet.pdf", gg1 + gg2,
        device = "pdf", height = 12, width = 12,
        path = paste0("figures/", age, "/byCellType/", cluster_name, "/cellchat/")
    )
    gg1c <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2))
    gg2c <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "weight")
    ggsave("compareInteractions.pdf", gg1c + gg2c,
        device = "pdf", height = 12, width = 12,
        path = paste0("figures/", age, "/byCellType/", cluster_name, "/cellchat/")
    )

    pdf(paste0("figures/", age, "/byCellType/", cluster_name, "/cellchat/newVisual_diffInteraction.pdf"), width = 8, heigh = 8)
    CellChat::netVisual_diffInteraction(cellchat, weight.scale = TRUE)
    dev.off()
    pdf(paste0("figures/", age, "/byCellType/", cluster_name, "/cellchat/newVisual_diffInteraction_weight.pdf"), width = 8, heigh = 8)
    CellChat::netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
    dev.off()
}

#' Get the best K value
#'
#' This function calculates the best K value for a plot by applying a weighted average
#' to the rate of change between consecutive values in two arrays. The weights are linearly increasing
#' and can be transformed by a factor defined by the user.
#'
#' @param plot A CellChat plot
#' @return The best K value
#' @export
GetBestK <- function(plot) {
    pg <- ggplot2::ggplot_build(plot)
    a <- pg$data[[1]][pg$data[[1]]$group == 1, ]$y
    b <- pg$data[[1]][pg$data[[1]]$group == 2, ]$y

    # Define the weights (linearly increasing weights)
    weights <- seq_along(a)

    transformation_function <- function(x, alpha) {
        return(x^alpha)
    }

    # Define the transformation factor
    transformation_factor <- 0.5

    # Normalize the arrays
    normalized_a <- (a - min(a)) / (max(a) - min(a))
    normalized_b <- (b - min(b)) / (max(b) - min(b))

    # Calculate the rate of change between consecutive values in each array
    diff_a <- c(0, diff(normalized_a))
    diff_b <- c(0, diff(normalized_b))

    # Apply the transformation function to the weights
    transformed_weights <- transformation_function(weights, transformation_factor)

    # Calculate the weighted average, considering the rate of change and transformed weights
    weighted_average <- (transformed_weights * normalized_a) + (transformed_weights * normalized_b)

    # Find the index with the maximum weighted average
    best_index <- which.max(weighted_average)
    return(best_index + 1)
}

SelectKParallel <- function(cellchat, pattern) {
    CellChat::selectK(cellchat, pattern = pattern)
}

#' Subset seurat object with less that n cells
#'
#'
#' @param seuratObj Seurat object with hyperoxia and Room air cells
#' @param cellNum minimum number of cells per cluster
#' @return Seurat object with less than n cells
#' @export
RemoveSmallClusters <- function(seuratObj, cellNum = 100) {
    rmCl <- table(seuratObj@active.ident) |>
        as.data.frame() |>
        tibble::column_to_rownames(var = "Var1") |>
        dplyr::filter(Freq < cellNum)
    seuratObj <- subset(seuratObj, ident = rownames(rmCl), invert = TRUE)
    return(seuratObj)
}
