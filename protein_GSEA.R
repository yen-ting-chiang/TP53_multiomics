################################################################################
# TP53 Protein-Level Preranked GSEA
#
# Uses moderated-t statistics from protein_differential_analysis/ as ranking
#
#   - Hallmark, C2 (CGP, BIOCARTA, KEGG_LEGACY, KEGG_MEDICUS, PID, REACTOME, WIKIPATHWAYS)
#   - C3 (TFT:GTRD, TFT:TFT_LEGACY), C4 (3CA, CGN, CM)
#   - C5 (GO:BP, GO:CC, GO:MF), C6 (Oncogenic), C9 (Cell Type)
#
# 9 Comparisons per dataset:
#   TP53mt vs TP53wt, MUT_GOF vs MUT_LOF, Hotspots vs MUT_LOF,
#   MUT_GOF vs TP53wt, MUT_LOF vs TP53wt, Hotspots vs TP53wt,
#   DN vs TP53wt, Non-DN vs TP53wt, DN vs non-DN
#
# Methodology: fgsea preranked GSEA with msigdbr gene sets
################################################################################

# ==============================================================================
# Section 1: Setup
# ==============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(msigdbr)
    library(fgsea)
    library(openxlsx)
})

set.seed(1234)

base_path <- "C:/Users/danny/Documents/R_project/TP53_multiomics"

datasets <- data.frame(
    folder = c("brca_cptac_2020", "coad_cptac_2019", "gbm_cptac_2021",
               "luad_cptac_2020", "paad_cptac_2021", "ucec_cptac_2020"),
    cancer_type = c("BRCA", "COAD", "GBM", "LUAD", "PAAD", "UCEC"),
    stringsAsFactors = FALSE
)

# Input: protein DEG results
deg_dir <- file.path(base_path, "protein_differential_analysis")

# Output: GSEA results
output_dir <- file.path(base_path, "protein_GSEA")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 9 comparisons with their file prefixes
comparisons <- data.frame(
    name = c("TP53mt_vs_TP53wt", "MUT_GOF_vs_MUT_LOF", "Hotspot_vs_MUT_LOF",
             "MUT_GOF_vs_TP53wt", "MUT_LOF_vs_TP53wt", "Hotspot_vs_TP53wt",
             "DN_vs_TP53wt", "NonDN_vs_TP53wt", "DN_vs_NonDN"),
    label = c("mt_vs_wt", "GOF_vs_LOF", "Hot_vs_LOF",
              "GOF_vs_wt", "LOF_vs_wt", "Hot_vs_wt",
              "DN_vs_wt", "NonDN_vs_wt", "DN_vs_NonDN"),
    stringsAsFactors = FALSE
)

# fgsea parameters
minSize <- 15
maxSize <- 500

cat("====================================================================\n")
cat("TP53 Protein-Level Preranked GSEA\n")
cat("====================================================================\n\n")

# ==============================================================================
# Section 2: Load Gene Sets from MSigDB
# ==============================================================================

cat("Loading MSigDB gene sets...\n")

# Helper: convert msigdbr output to named list for fgsea
msigdbr_to_list <- function(df) {
    split(df$gene_symbol, df$gs_name)
}

collections_to_fetch <- list(
    list(name="Hallmark",               cat="H",  sub=NULL),
    list(name="C2_CGP",                 cat="C2", sub="CGP"),
    list(name="C2_CP_BIOCARTA",         cat="C2", sub="CP:BIOCARTA"),
    list(name="C2_CP_KEGG_LEGACY",      cat="C2", sub="CP:KEGG_LEGACY"),
    list(name="C2_CP_KEGG_MEDICUS",     cat="C2", sub="CP:KEGG_MEDICUS"),
    list(name="C2_CP_PID",              cat="C2", sub="CP:PID"),
    list(name="C2_CP_REACTOME",         cat="C2", sub="CP:REACTOME"),
    list(name="C2_CP_WIKIPATHWAYS",     cat="C2", sub="CP:WIKIPATHWAYS"),
    list(name="C3_TFT_GTRD",            cat="C3", sub="TFT:GTRD"),
    list(name="C3_TFT_TFT_LEGACY",      cat="C3", sub="TFT:TFT_LEGACY"),
    list(name="C4_3CA",                 cat="C4", sub="3CA"),
    list(name="C4_CGN",                 cat="C4", sub="CGN"),
    list(name="C4_CM",                  cat="C4", sub="CM"),
    list(name="C5_GO_BP",               cat="C5", sub="GO:BP"),
    list(name="C5_GO_CC",               cat="C5", sub="GO:CC"),
    list(name="C5_GO_MF",               cat="C5", sub="GO:MF"),
    list(name="C6_Oncogenic",           cat="C6", sub=NULL),
    list(name="C9_CellType",            cat="C9", sub=NULL)
)

collections <- list()

for (info in collections_to_fetch) {
    if (is.null(info$sub)) {
        df <- suppressMessages(msigdbr(species = "Homo sapiens", category = info$cat))
    } else {
        df <- suppressMessages(msigdbr(species = "Homo sapiens", category = info$cat, subcategory = info$sub))
    }
    df <- df %>% dplyr::select(gs_name, gene_symbol)
    collections[[info$name]] <- msigdbr_to_list(df)
    cat(sprintf("  %-22s : %d gene sets\n", info$name, length(collections[[info$name]])))
}
cat("\n")

# ==============================================================================
# Section 3: Run Preranked GSEA
# ==============================================================================

#' Run fgsea preranked GSEA from a DEG table
#' @param deg_file Path to DEG CSV file
#' @param pathways Named list of gene sets
#' @return fgsea result data.frame or NULL
run_preranked_gsea <- function(deg_file, pathways) {
    if (!file.exists(deg_file)) return(NULL)

    deg <- fread(deg_file)
    if (!("gene" %in% colnames(deg)) || !("t" %in% colnames(deg))) {
        cat("    [WARN] Missing gene/t columns in:", basename(deg_file), "\n")
        return(NULL)
    }

    # Build named vector of moderated-t statistics
    stats_vec <- setNames(deg$t, deg$gene)

    # Remove duplicates (keep first = most significant from topTable)
    stats_vec <- stats_vec[!duplicated(names(stats_vec))]

    # Remove non-finite values
    stats_vec <- stats_vec[is.finite(stats_vec)]

    if (length(stats_vec) < 100) {
        cat("    [WARN] Too few genes:", length(stats_vec), "\n")
        return(NULL)
    }

    # Sort descending (required by fgsea)
    stats_vec <- sort(stats_vec, decreasing = TRUE)

    # Run fgsea
    res <- tryCatch({
        suppressWarnings(fgsea::fgseaMultilevel(
            pathways = pathways,
            stats = stats_vec,
            minSize = minSize,
            maxSize = maxSize,
            eps = 0
        ))
    }, error = function(e) {
        cat("    [ERROR] fgsea failed:", conditionMessage(e), "\n")
        NULL
    })

    if (is.null(res) || nrow(res) == 0) return(NULL)

    # Convert leadingEdge list to semicolon-separated string for CSV/XLSX output
    res$leadingEdge <- vapply(res$leadingEdge,
                              function(x) paste(x, collapse = ";"),
                              character(1))

    # Sort by padj
    res <- res[order(res$padj), ]
    as.data.frame(res)
}

# ==============================================================================
# Section 4: Main Processing Loop
# ==============================================================================

# Storage for summary statistics (one per collection)
all_summaries <- list()
for (coll_name in names(collections)) {
    all_summaries[[coll_name]] <- list()
}

for (i in seq_len(nrow(datasets))) {
    ds_folder <- datasets$folder[i]
    ds_cancer <- datasets$cancer_type[i]

    cat("==================================================================\n")
    cat("Processing:", ds_folder, "(", ds_cancer, ")\n")
    cat("==================================================================\n")

    ds_deg_dir <- file.path(deg_dir, ds_folder)
    if (!dir.exists(ds_deg_dir)) {
        cat("  [SKIP] DEG directory not found\n\n")
        next
    }

    ds_out <- file.path(output_dir, ds_folder)
    dir.create(ds_out, recursive = TRUE, showWarnings = FALSE)

    for (coll_name in names(collections)) {
        pathways <- collections[[coll_name]]
        cat("\n  --- Collection:", coll_name, "---\n")

        # Create collection output subdirectory
        coll_out <- file.path(ds_out, coll_name)
        dir.create(coll_out, recursive = TRUE, showWarnings = FALSE)

        # Summary row for this dataset + collection
        summary_row <- list(dataset = ds_folder, cancer_type = ds_cancer)

        for (j in seq_len(nrow(comparisons))) {
            comp_name <- comparisons$name[j]
            comp_label <- comparisons$label[j]

            deg_file <- file.path(ds_deg_dir, paste0("DEG_", comp_name, ".csv"))
            cat("    ", comp_label, ":", sep = "")

            res <- run_preranked_gsea(deg_file, pathways)

            if (!is.null(res)) {
                n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
                n_up <- sum(res$padj < 0.05 & res$NES > 0, na.rm = TRUE)
                n_down <- sum(res$padj < 0.05 & res$NES < 0, na.rm = TRUE)
                cat(" ", nrow(res), "sets tested,", n_sig, "significant\n")

                # Save CSV
                fwrite(res, file.path(coll_out, paste0("GSEA_", comp_name, ".csv")))

                # Save XLSX
                wb <- createWorkbook()
                addWorksheet(wb, comp_label)
                writeData(wb, 1, res)
                saveWorkbook(wb, file.path(coll_out, paste0("GSEA_", comp_name, ".xlsx")),
                             overwrite = TRUE)

                summary_row[[paste0(comp_name, "_total")]] <- nrow(res)
                summary_row[[paste0(comp_name, "_sig")]] <- n_sig
                summary_row[[paste0(comp_name, "_up")]] <- n_up
                summary_row[[paste0(comp_name, "_down")]] <- n_down
            } else {
                cat(" no DEG file or insufficient data\n")
                summary_row[[paste0(comp_name, "_total")]] <- NA
                summary_row[[paste0(comp_name, "_sig")]] <- NA
                summary_row[[paste0(comp_name, "_up")]] <- NA
                summary_row[[paste0(comp_name, "_down")]] <- NA
            }
        }

        all_summaries[[coll_name]][[ds_folder]] <-
            as.data.frame(summary_row, stringsAsFactors = FALSE)
    }
    cat("\n")
}

# ==============================================================================
# Section 5: Summary Statistics (per collection)
# ==============================================================================

cat("====================================================================\n")
cat("Generating Summary Statistics\n")
cat("====================================================================\n\n")

for (coll_name in names(all_summaries)) {
    cat("--- Summary for", coll_name, "---\n")

    summary_df <- bind_rows(all_summaries[[coll_name]])

    # CSV
    fwrite(summary_df, file.path(output_dir,
        paste0("GSEA_summary_", coll_name, ".csv")))

    # XLSX
    wb_sum <- createWorkbook()
    addWorksheet(wb_sum, "Summary")
    writeData(wb_sum, "Summary", summary_df)
    header_style <- createStyle(
        textDecoration = "bold", halign = "center",
        border = "bottom", fgFill = "#4472C4", fontColour = "white"
    )
    addStyle(wb_sum, "Summary", header_style,
             rows = 1, cols = 1:ncol(summary_df), gridExpand = TRUE)
    setColWidths(wb_sum, "Summary", cols = 1:ncol(summary_df), widths = "auto")
    saveWorkbook(wb_sum, file.path(output_dir,
        paste0("GSEA_summary_", coll_name, ".xlsx")), overwrite = TRUE)

    # Print significant counts only
    sig_cols <- grep("_sig$", colnames(summary_df), value = TRUE)
    print_df <- summary_df[, c("dataset", "cancer_type", sig_cols)]
    print(as.data.frame(print_df), row.names = FALSE)
    cat("\n")
}

cat("====================================================================\n")
cat("Pipeline Complete!\n")
cat("All results saved to:", output_dir, "\n")
cat("====================================================================\n")
