#' Calculate enrichment scores for a given gene set and group
#'
#' This function takes a gene set and three data frames containing gene set enrichment analysis results for a given group and calculates
#' enrichment scores for each term in the gene set. The function first subsets the data frames by terms with a significant enrichment (p-value < 0.05)
#' and calculates the negative log10 of the p-value. The function then subsets the data frames by terms with at least 5 overlapping genes and counts
#' the number of overlapping genes. The function returns a list containing three data frames with the enrichment scores for the upregulated, downregulated,
#' and all genes in the group.
#'
#' @param name A character string specifying the name of the group.
#' @param oxyUp A data frame containing gene set enrichment analysis results for the upregulated genes in the group.
#' @param oxyDown A data frame containing gene set enrichment analysis results for the downregulated genes in the group.
#' @param oxyAll A data frame containing gene set enrichment analysis results for all genes in the group.
#'
#' @return A list containing three data frames with the enrichment scores for the upregulated, downregulated, and all genes in the group.
#'
#' @export
GetScores <- function(name, oxyUp, oxyDown, oxyAll) {
    geneSetUp <- oxyUp %>%
        dplyr::filter(P.value < 0.05) %>%
        dplyr::mutate(log10Pval = round(-log10(P.value), 3)) %>%
        dplyr::select(Term, log10Pval, Genes) %>%
        dplyr::mutate(Gene.Count = sapply(
            stringr::str_split(Genes, ";"),
            function(x) length(unique(x))
        )) %>%
        dplyr::filter(Gene.Count >= 5) %>%
        tibble::column_to_rownames("Term") %>%
        setNames(c(
            paste0("-log10(pVal) ", name, "_UP"),
            paste0("Overlapping genes ", name, "_UP"),
            paste0("Overlapping gene number ", name, "_UP")
        ))
    geneSetDown <- oxyDown %>%
        dplyr::filter(P.value < 0.05) %>%
        dplyr::mutate(log10Pval = round(-log10(P.value), 3)) %>%
        dplyr::select(Term, log10Pval, Genes) %>%
        dplyr::mutate(Gene.Count = sapply(
            stringr::str_split(Genes, ";"),
            function(x) length(unique(x))
        )) %>%
        dplyr::filter(Gene.Count >= 5) %>%
        tibble::column_to_rownames("Term") %>%
        setNames(c(
            paste0("-log10(pVal) ", name, "_DOWN"),
            paste0("Overlapping genes ", name, "_DOWN"),
            paste0("Overlapping gene number ", name, "_DOWN")
        ))
    geneSetAll <- oxyAll %>%
        dplyr::filter(P.value < 0.05) %>%
        dplyr::mutate(log10Pval = round(-log10(P.value), 3)) %>%
        dplyr::select(Term, log10Pval, Genes) %>%
        dplyr::mutate(Gene.Count = sapply(
            stringr::str_split(Genes, ";"),
            function(x) length(unique(x))
        )) %>%
        dplyr::filter(Gene.Count >= 5) %>%
        tibble::column_to_rownames("Term") %>%
        setNames(c(
            paste0("-log10(pVal) ", name, "_ALL"),
            paste0("Overlapping genes ", name, "_ALL"),
            paste0("Overlapping gene number ", name, "_ALL")
        ))
    return(list(geneSetUp, geneSetDown, geneSetAll))
}

#' Calculate enrichment scores for a given gene set and group using GSEA
#'
#' This function takes a gene set and three data frames containing gene set enrichment analysis results for a given group and calculates
#' enrichment scores for each term in the gene set using the GSEA algorithm. The function first subsets the data frames by terms with a
#' significant enrichment (p-value < 0.05) and calculates the normalized enrichment score (NES). The function then subsets the data frames
#' by terms with at least 5 overlapping genes and counts the number of overlapping genes. The function returns a list containing three data
#' frames with the enrichment scores for the upregulated, downregulated, and all genes in the group.
#'
#' @param name A character string specifying the name of the group.
#' @param oxyUp A data frame containing gene set enrichment analysis results for the upregulated genes in the group.
#' @param oxyDown A data frame containing gene set enrichment analysis results for the downregulated genes in the group.
#' @param oxyAll A data frame containing gene set enrichment analysis results for all genes in the group.
#' @param sets A named list of gene sets to be used in the GSEA analysis.
#'
#' @return A list containing three data frames with the enrichment scores for the upregulated, downregulated, and all genes in the group.
#'
#' @export
GetScoresGSEA <- function(name, oxyUp, oxyDown, oxyAll, sets) {
    fgseaResUp <- tryCatch(
        {
            fgsea::fgsea(sets, stats = tibble::deframe(oxyUp)) %>%
                tibble::as_tibble() %>%
                dplyr::filter(pval < 0.05) %>%
                dplyr::mutate(NES = round(NES, 3)) %>%
                dplyr::arrange(desc(NES)) %>%
                dplyr::select(pathway, NES, size, leadingEdge) %>%
                dplyr::filter(size >= 5) %>%
                tibble::column_to_rownames("pathway") %>%
                setNames(c(paste0("NES ", name, "_UP"), paste0("Overlapping gene number ", name, "_UP"), paste0("Overlapping genes ", name, "_UP")))
        },
        error = function(e) {
            data.frame(matrix(ncol = 3, nrow = 0, dimnames = list(NULL, c(
                paste0("NES ", name, "_UP"),
                paste0("Overlapping gene number ", name, "_UP"),
                paste0("Overlapping genes ", name, "_UP")
            ))), check.names = FALSE)
        }
    )
    fgseaResDown <- tryCatch(
        {
            fgsea::fgsea(sets, stats = tibble::deframe(oxyDown)) %>%
                tibble::as_tibble() %>%
                dplyr::filter(pval < 0.05) %>%
                dplyr::mutate(NES = round(NES, 3)) %>%
                dplyr::arrange(desc(NES)) %>%
                dplyr::select(pathway, NES, size, leadingEdge) %>%
                dplyr::filter(size >= 5) %>%
                tibble::column_to_rownames("pathway") %>%
                setNames(c(paste0("NES ", name, "_DOWN"), paste0("Overlapping gene number ", name, "_DOWN"), paste0("Overlapping genes ", name, "_DOWN")))
        },
        error = function(e) {
            data.frame(matrix(ncol = 3, nrow = 0, dimnames = list(NULL, c(
                paste0("NES ", name, "_DOWN"),
                paste0("Overlapping gene number ", name, "_DOWN"),
                paste0("Overlapping genes ", name, "_DOWN")
            ))), check.names = FALSE)
        }
    )
    fgseaResAll <- tryCatch(
        {
            fgsea::fgsea(sets, stats = tibble::deframe(oxyAll)) %>%
                tibble::as_tibble() %>%
                dplyr::filter(pval < 0.05) %>%
                dplyr::mutate(NES = round(NES, 3)) %>%
                dplyr::arrange(desc(NES)) %>%
                dplyr::select(pathway, NES, size, leadingEdge) %>%
                dplyr::filter(size >= 5) %>%
                tibble::column_to_rownames("pathway") %>%
                setNames(c(
                    paste0("NES ", name, "_ALL"),
                    paste0("Overlapping gene number ", name, "_ALL"),
                    paste0("Overlapping genes ", name, "_ALL")
                ))
        },
        error = function(e) {
            data.frame(matrix(ncol = 3, nrow = 0, dimnames = list(NULL, c(
                paste0("NES ", name, "_ALL"),
                paste0("Overlapping gene number ", name, "_ALL"),
                paste0("Overlapping genes ", name, "_ALL")
            ))), check.names = FALSE)
        }
    )
    return(list(fgseaResUp, fgseaResDown, fgseaResAll))
}

#' Generate a final matrix of gene set enrichment results
#'
#' This function takes a gene set matrix and generates a final matrix of gene set enrichment results.
#' The function first selects the gene sets from the gene set matrix and subsets them into three data
#' frames containing the enrichment scores, the number of overlapping genes, and the overlapping gene names.
#' The function then joins the three data frames into a single data frame and returns it.
#'
#' @param gene_matrix A data frame containing gene set enrichment results.
#' @param start_string A character string specifying the prefix of the column names to be selected.
#'
#' @return A data frame containing the gene set enrichment results.
#'
#' @export
GetFinalMatrix <- function(gene_matrix, start_string = "-log10") {
    dbGenes <- gene_matrix %>%
        purrr::map(~ .x %>%
            tibble::rownames_to_column("rn")) %>%
        purrr::reduce(dplyr::full_join, by = "rn") %>%
        tibble::column_to_rownames("rn")
    db_log10 <- dbGenes %>% dplyr::select(starts_with(start_string))
    db_gen_num <- dbGenes %>% dplyr::select(starts_with("Overlapping gene num"))
    db_gen_name <- dbGenes %>% dplyr::select(starts_with("Overlapping genes"))

    dbGenes <- list(db_log10, db_gen_num, db_gen_name) %>%
        purrr::map(~ .x %>% tibble::rownames_to_column("rn")) %>%
        purrr::reduce(dplyr::full_join, by = "rn") %>%
        tibble::column_to_rownames("rn") %>%
        tibble::rownames_to_column("GeneSet")

    if (nrow(dbGenes) == 0) {
        dbGenes[1, ] <- rep(0, ncol(dbGenes))
    }

    return(dbGenes)
}

#' Generate enrichment summary for a single condition
#'
#' This function takes a list of differentially expressed genes for each cluster, as well as information about the cell type and age,
#' and generates an enrichment summary for a single condition. The function uses the Enrichr API to perform gene set enrichment analysis
#' and generates a summary of the results. The summary is saved to a file with a name that includes the cluster name, cell type, age,
#' and an optional string to append to the file name.
#'
#' @param clusterDE A list of differentially expressed genes for each cluster.
#' @param cellType A character string specifying the cell type.
#' @param age A character string specifying the age.
#' @param appendToFileName A character string to append to the file name.
#' @param condition A character string specifying the condition.
#'
#' @return None.
#'
#' @export
GetEnrichmentSummarySingleCondition <- function(clusterDE, cellType, age, appendToFileName = "", condition = "bySex") {
    options(enrichR.base.address = "https://maayanlab.cloud/Enrichr/")
    mDf <- msigdbr::msigdbr(species = "Mus musculus", category = "C5")
    fgseaSets <- mDf %>% split(x = .$gene_symbol, f = .$gs_name)
    for (clusterName in names(clusterDE)) {
        DEGenes <- clusterDE[[clusterName]]
        tryCatch(
            {
                GenerateSummarySingleCondition(clusterName, DEGenes, cellType, age, fgseaSets, appendToFileName, condition)
            },
            error = function(e) {
                cat("ERROR :", conditionMessage(e), "\n")
            }
        )
    }
}

#' Generate enrichment summary for a single condition and cluster
#'
#' This function takes a list of differentially expressed genes for a single cluster, as well as information about the cell type and age,
#' and generates an enrichment summary for the cluster. The function uses the Enrichr API to perform gene set enrichment analysis and
#' generates a summary of the results. The summary is saved to a file with a name that includes the cluster name, cell type, age, and an
#' optional string to append to the file name.
#'
#' @param clusterName A character string specifying the name of the cluster.
#' @param DEGenes A character vector containing the differentially expressed genes for the cluster.
#' @param cellType A character string specifying the cell type.
#' @param age A character string specifying the age.
#' @param fgseaSets A list of gene sets to use for enrichment analysis.
#' @param appendToFileName A character string to append to the file name.
#' @param condition A character string specifying the condition.
#'
#' @return None.
#'
#' @export
GenerateSummarySingleCondition <- function(clusterName, DEGenes, cellType, age, fgseaSets, appendToFileName = "", condition = "bySex") {
    # Spread sheet styling
    colStyle <- openxlsx::createStyle(textRotation = 45, textDecoration = "bold")
    rBorder <- openxlsx::createStyle(border = "Right")
    bBorder <- openxlsx::createStyle(border = "Bottom")
    notEmpty <- openxlsx::createStyle(bgFill = "#e8e8e8")
    boldText <- openxlsx::createStyle(textDecoration = "bold")
    headerText <- openxlsx::createStyle(fontSize = 18, textDecoration = "bold")
    wb <- openxlsx::createWorkbook()

    # load all datasets needed
    maleVsFemale <- DEGenes %>% dplyr::filter(sig_exp != "NO")

    # organize gene sets into a list.
    sections <- list(
        list("MvsF", maleVsFemale)
    )

    # list of analysis to run
    nameSets <- c("Hallmark Combined", "Kegg Combined", "Reactome Combined", "GOBP Combined", "GSEA Combined")
    # intializize empty matrix
    geneNum <- as.data.frame(matrix(nrow = length(sections) * 3, ncol = length(nameSets) + 1))
    # organize the column main page
    colnames(geneNum) <- c("# DEGs", "Hallmark", "Kegg", "Reactome", "GOBP", "GSEA")
    geneSets <- c(
        "MvsF"
    )
    setsColors <- c("#f27f0c", "#de43f0", "#4d7cff", "#de43f0", "#4d7cff", "#de43f0", "#4d7cff")

    # generate gene set names
    geneSetsNames <- c()
    for (name in geneSets) {
        geneSetsNames <- c(geneSetsNames, paste0("-log10(pVal) ", name, "_UP"))
        geneSetsNames <- c(geneSetsNames, paste0("-log10(pVal) ", name, "_DOWN"))
        geneSetsNames <- c(geneSetsNames, paste0("-log10(pVal) ", name, "_ALL"))
    }
    rownames(geneNum) <- geneSetsNames

    hallmarkGenes <- list()
    keggGenes <- list()
    reactomeGenes <- list()
    GOBPGenes <- list()
    GSEAGenes <- list()
    # for every gene set get the enrichment sets
    for (i in sections) {
        dbs <- c(
            "MSigDB_Hallmark_2020",
            "KEGG_2019_Mouse",
            "Reactome_2016",
            "GO_Biological_Process_2021"
        )

        # empty pathway list in case no genes were given
        edf <- data.frame(matrix(ncol = 9, nrow = 0, dimnames = list(NULL, c(
            "Term", "Overlap", "P.value",
            "Adjusted.P.value", "Old.P.value", "Old.Adjusted.P.value",
            "Odds.Ratio", "Combined.Score", "Genes"
        ))))
        epth <- list(
            MSigDB_Hallmark_2020 = edf, KEGG_2019_Mouse = edf,
            Reactome_2016 = edf, GO_Biological_Process_2021 = edf
        )

        markers <- i[[2]]
        set <- rownames(markers[markers$sig_exp == "UP", ])
        geneNum[paste0("-log10(pVal) ", unlist(i[[1]]), "_UP"), "# DEGs"] <- length(set)
        capture.output(markersUp <- enrichR::enrichr(set, dbs))
        markersUp <- if (is.null(markersUp)) epth else markersUp
        set <- rownames(markers[markers$sig_exp == "DOWN", ])
        geneNum[paste0("-log10(pVal) ", unlist(i[[1]]), "_DOWN"), "# DEGs"] <- length(set)
        capture.output(markersDown <- enrichR::enrichr(set, dbs))
        markersDown <- if (is.null(markersDown)) epth else markersDown
        set <- c(rownames(markers[markers$sig_exp == "UP", ]), rownames(markers[markers$sig_exp == "DOWN", ]))
        geneNum[paste0("-log10(pVal) ", unlist(i[[1]]), "_ALL"), "# DEGs"] <- length(set)
        capture.output(markersAll <- enrichR::enrichr(set, dbs))
        markersAll <- if (is.null(markersAll)) epth else markersAll

        hallmarkGenes <- c(hallmarkGenes, GetScores(
            unlist(i[[1]]), markersUp$MSigDB_Hallmark_2020,
            markersDown$MSigDB_Hallmark_2020,
            markersAll$MSigDB_Hallmark_2020
        ))
        keggGenes <- c(keggGenes, GetScores(
            unlist(i[[1]]), markersUp$KEGG_2019_Mouse,
            markersDown$KEGG_2019_Mouse,
            markersAll$KEGG_2019_Mouse
        ))
        reactomeGenes <- c(reactomeGenes, GetScores(
            unlist(i[[1]]), markersUp$Reactome_2016,
            markersDown$Reactome_2016,
            markersAll$Reactome_2016
        ))
        GOBPGenes <- c(GOBPGenes, GetScores(
            unlist(i[[1]]), markersUp$GO_Biological_Process_2021,
            markersDown$GO_Biological_Process_2021,
            markersAll$GO_Biological_Process_2021
        ))

        gseaUp <- markers[markers$sig_exp == "UP", ] %>%
            tibble::rownames_to_column("feature") %>%
            dplyr::arrange(desc(avg_log2FC)) %>%
            dplyr::select(feature, avg_log2FC)
        gseaDown <- markers[markers$sig_exp == "DOWN", ] %>%
            tibble::rownames_to_column("feature") %>%
            dplyr::arrange(desc(avg_log2FC)) %>%
            dplyr::select(feature, avg_log2FC)
        gseaAll <- markers[markers$sig_exp == "UP" | markers$sig_exp == "DOWN", ] %>%
            tibble::rownames_to_column("feature") %>%
            dplyr::arrange(desc(avg_log2FC)) %>%
            dplyr::select(feature, avg_log2FC)

        GSEAGenes <- c(GSEAGenes, GetScoresGSEA(unlist(i[[1]]), gseaUp, gseaDown, gseaAll, fgseaSets))
    }

    # get final data frame to write to excel
    hallmarkGenes <- GetFinalMatrix(hallmarkGenes)
    keggGenes <- GetFinalMatrix(keggGenes)
    reactomeGenes <- GetFinalMatrix(reactomeGenes)
    GOBPGenes <- GetFinalMatrix(GOBPGenes)
    GSEAGenes <- GetFinalMatrix(GSEAGenes, start_string = "NES")

    listH <- list(nrow(hallmarkGenes), nrow(keggGenes), nrow(reactomeGenes), nrow(GOBPGenes), nrow(GSEAGenes))
    listW <- list(ncol(hallmarkGenes), ncol(keggGenes), ncol(reactomeGenes), ncol(GOBPGenes), ncol(GSEAGenes))
    names(listH) <- nameSets
    names(listW) <- nameSets

    # add all needed worksheets and write data frames in them
    openxlsx::addWorksheet(wb, sheetName = "summary", gridLines = TRUE, zoom = 120)
    for (set in nameSets) {
        openxlsx::addWorksheet(wb, sheetName = set, gridLines = TRUE, zoom = 100)
        openxlsx::freezePane(wb, sheet = set, firstRow = TRUE, firstCol = TRUE)
    }
    openxlsx::writeData(wb, x = hallmarkGenes, sheet = "Hallmark Combined", startRow = 1, rowNames = FALSE)
    openxlsx::writeData(wb, x = keggGenes, sheet = "Kegg Combined", startRow = 1, rowNames = FALSE)
    openxlsx::writeData(wb, x = reactomeGenes, sheet = "Reactome Combined", startRow = 1, rowNames = FALSE)
    openxlsx::writeData(wb, x = GOBPGenes, sheet = "GOBP Combined", startRow = 1, rowNames = FALSE)
    openxlsx::writeData(wb, x = GSEAGenes, sheet = "GSEA Combined", startRow = 1, rowNames = FALSE)

    # add style to every gene set
    for (set in nameSets) {
        height <- listH[[set]]
        width <- listW[[set]]
        openxlsx::addStyle(wb, set,
            cols = 1:width, rows = 1, style = colStyle,
            gridExpand = TRUE, stack = TRUE
        )
        openxlsx::addStyle(wb, set,
            cols = 1:width, rows = 1, style = bBorder,
            gridExpand = TRUE, stack = TRUE
        )
        openxlsx::addStyle(wb, set,
            cols = c(1, seq(4, width, by = 3)), rows = 1:height,
            style = rBorder, gridExpand = TRUE, stack = TRUE
        )
        openxlsx::conditionalFormatting(wb, set,
            cols = 2:(length(sections) * 3),
            rows = 2:height,
            style = c("#ffffff", "#ff9999"),
            type = "colourScale"
        )
        openxlsx::conditionalFormatting(wb, set,
            cols = ((length(sections) * 6) + 2):((length(sections) * 3) * 3 + 1),
            rows = 2:height, rule = "!=0", style = notEmpty
        )
        openxlsx::setColWidths(wb, set, cols = 1, widths = 35)
        i <- 2
        for (color in c(setsColors, setsColors, setsColors)) {
            openxlsx::addStyle(wb, set,
                cols = i:(i + 2), rows = 1,
                style = openxlsx::createStyle(fontColour = color),
                gridExpand = TRUE, stack = TRUE
            )
            i <- i + 3
        }
    }

    # get total number of genes per column
    totalGenes <- as.data.frame(matrix(nrow = 1, ncol = length(nameSets)))
    i <- 2
    for (geneSet in list(hallmarkGenes, keggGenes, reactomeGenes, GOBPGenes)) {
        for (colname in rownames(geneNum)) {
            geneNum[colname, i] <- length(na.omit(geneSet[, colname]))
        }
        totalGenes[1, i - 1] <- nrow(geneSet)
        i <- i + 1
    }
    for (colname in rownames(geneNum)) {
        cl <- gsub("-log10(pVal)", "NES", colname, fixed = TRUE)
        geneNum[colname, i] <- length(na.omit(GSEAGenes[, cl]))
    }
    totalGenes[1, i - 1] <- nrow(GSEAGenes)

    openxlsx::writeData(wb, x = geneNum, sheet = "summary", startRow = 3, startCol = 2, rowNames = TRUE)
    openxlsx::writeData(wb, x = totalGenes, sheet = "summary", startRow = 2, startCol = 4, rowNames = FALSE, colNames = FALSE)
    openxlsx::writeData(wb, "summary", "Number of pathways enriched", startCol = 4, startRow = 1)
    openxlsx::mergeCells(wb, "summary", cols = 4:8, rows = 1)

    openxlsx::addStyle(wb, "summary",
        cols = 2:(2 + ncol(geneNum)), rows = seq(3, 3 + nrow(geneNum), by = 3), style = bBorder,
        gridExpand = TRUE, stack = TRUE
    )
    openxlsx::addStyle(wb, "summary",
        cols = 3:(2 + ncol(geneNum)), rows = 2, style = bBorder,
        gridExpand = TRUE, stack = TRUE
    )
    openxlsx::addStyle(wb, "summary",
        cols = c(2, 3, 8), rows = seq(3, 3 + nrow(geneNum)), style = rBorder,
        gridExpand = TRUE, stack = TRUE
    )
    openxlsx::addStyle(wb, "summary",
        cols = 1, rows = seq(4, 3 + nrow(geneNum)), style = rBorder,
        gridExpand = TRUE, stack = TRUE
    )
    openxlsx::addStyle(wb, "summary",
        cols = 4:(2 + ncol(geneNum)), rows = 3, style = boldText,
        gridExpand = TRUE, stack = TRUE
    )
    openxlsx::addStyle(wb, "summary",
        cols = 4:(2 + ncol(geneNum)), rows = 1, style = headerText,
        gridExpand = TRUE, stack = TRUE
    )
    openxlsx::conditionalFormatting(wb, "summary",
        cols = 4:(3 + length(nameSets)),
        rows = 4:(3 + length(sections) * 3),
        style = c("#ffffff", "#7be380"),
        type = "colourScale"
    )
    openxlsx::setColWidths(wb, "summary", cols = 2, widths = "auto")

    i <- 4
    for (color in c(setsColors, setsColors, setsColors)) {
        openxlsx::addStyle(wb, "summary",
            cols = 2, rows = i:(i + 2),
            style = openxlsx::createStyle(fontColour = color),
            gridExpand = TRUE, stack = TRUE
        )
        i <- i + 3
    }
    openxlsx::saveWorkbook(wb, paste0(
        "./results/", age, "/byCellType/", cellType, "/",
        condition, "/", clusterName, "/", clusterName, "_enrichment_summary", appendToFileName, ".xlsx"
    ), overwrite = TRUE)
}
