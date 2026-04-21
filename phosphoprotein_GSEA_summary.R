################################################################################
# TP53 Phosphoprotein GSEA Results Summarization
#
# Aggregates the results from the fgsea analysis (phosphoprotein_GSEA)
#
# For each dataset and collection across the 9 comparisons, this script tallies for each pathway:
#   - Total number of times it is significant (padj < 0.05)
#   - Total positive significant (padj < 0.05 & NES > 0)
#   - Total negative significant (padj < 0.05 & NES < 0)
#   - A list of the specific comparisons where it was positive significant.
#   - A list of the specific comparisons where it was negative significant.
#
# Outputs: .csv and .xlsx files saved into phosphoprotein_GSEA folder.
################################################################################

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(openxlsx)
})

base_path <- "C:/Users/danny/Documents/R_project/TP53_multiomics"
gsea_dir <- file.path(base_path, "phosphoprotein_GSEA")

datasets <- data.frame(
    folder = c("brca_cptac_2020", "coad_cptac_2019", "gbm_cptac_2021",
               "luad_cptac_2020", "paad_cptac_2021", "ucec_cptac_2020"),
    cancer_type = c("BRCA", "COAD", "GBM", "LUAD", "PAAD", "UCEC"),
    stringsAsFactors = FALSE
)

comparisons <- data.frame(
    name = c("TP53mt_vs_TP53wt", "MUT_GOF_vs_MUT_LOF", "Hotspot_vs_MUT_LOF",
             "MUT_GOF_vs_TP53wt", "MUT_LOF_vs_TP53wt", "Hotspot_vs_TP53wt",
             "DN_vs_TP53wt", "NonDN_vs_TP53wt", "DN_vs_NonDN"),
    label = c("mt_vs_wt", "GOF_vs_LOF", "Hot_vs_LOF",
              "GOF_vs_wt", "LOF_vs_wt", "Hot_vs_wt",
              "DN_vs_wt", "NonDN_vs_wt", "DN_vs_NonDN"),
    stringsAsFactors = FALSE
)

collections <- c("PTMsigDB")

cat("====================================================================\n")
cat("TP53 Phosphoprotein GSEA Summarization\n")
cat("====================================================================\n\n")

if (!dir.exists(gsea_dir)) {
    stop("GSEA directory not found. Please run the GSEA script first.")
}

for (i in seq_len(nrow(datasets))) {
    ds_folder <- datasets$folder[i]
    ds_cancer <- datasets$cancer_type[i]

    cat("Processing Summary for Dataset:", ds_folder, "(", ds_cancer, ")\n")

    ds_out <- file.path(gsea_dir, ds_folder)
    if (!dir.exists(ds_out)) {
        cat("  [SKIP] Dataset GSEA directory not found\n")
        next
    }

    for (coll in collections) {
        coll_dir <- file.path(ds_out, coll)
        if (!dir.exists(coll_dir)) {
            next
        }

        # Gather data across all 9 comparisons for this collection
        all_pathways <- list()

        for (j in seq_len(nrow(comparisons))) {
            comp_name <- comparisons$name[j]
            comp_label <- comparisons$label[j]

            gsea_file <- file.path(coll_dir, paste0("GSEA_", comp_name, ".csv"))

            if (file.exists(gsea_file)) {
                gsea_res <- tryCatch(fread(gsea_file), error = function(e) NULL)
                
                if (!is.null(gsea_res) && nrow(gsea_res) > 0 && "pathway" %in% colnames(gsea_res)) {
                    # Format specific df for binding
                    df <- data.frame(
                        pathway = gsea_res$pathway,
                        padj = gsea_res$padj,
                        NES = gsea_res$NES,
                        comparison = comp_label,
                        stringsAsFactors = FALSE
                    )
                    all_pathways[[comp_name]] <- df
                }
            }
        }

        if (length(all_pathways) == 0) {
            cat("  [SKIP] No GSEA results found for collection:", coll, "\n")
            next
        }

        # Bind all rows
        combined_df <- bind_rows(all_pathways)

        # Summarize per pathway
        summary_df <- combined_df %>%
            group_by(pathway) %>%
            summarize(
                Total_Tested = n(),
                Sig_Count = sum(!is.na(padj) & padj < 0.05),
                Pos_Sig_Count = sum(!is.na(padj) & padj < 0.05 & !is.na(NES) & NES > 0),
                Neg_Sig_Count = sum(!is.na(padj) & padj < 0.05 & !is.na(NES) & NES < 0),
                Pos_Significant_In = paste(comparison[!is.na(padj) & padj < 0.05 & !is.na(NES) & NES > 0], collapse = ", "),
                Neg_Significant_In = paste(comparison[!is.na(padj) & padj < 0.05 & !is.na(NES) & NES < 0], collapse = ", ")
            ) %>%
            arrange(desc(Sig_Count), pathway)

        # Filter to only keep pathways that are strictly significant at least once
        summary_df_sig <- summary_df %>% filter(Sig_Count > 0)

        # Write per-dataset outputs for this collection
        csv_out <- file.path(coll_dir, paste0("Detailed_Pathway_Summary_", ds_folder, "_", coll, ".csv"))
        fwrite(summary_df, csv_out)
        
        xlsx_out <- file.path(coll_dir, paste0("Detailed_Pathway_Summary_", ds_folder, "_", coll, ".xlsx"))
        wb <- createWorkbook()
        
        addWorksheet(wb, "All_Tested_Pathways")
        writeData(wb, "All_Tested_Pathways", summary_df)
        
        addWorksheet(wb, "Significant_Only")
        writeData(wb, "Significant_Only", summary_df_sig)
        
        # Styles for Excel
        header_style <- createStyle(
            textDecoration = "bold", halign = "center",
            border = "bottom", fgFill = "#4472C4", fontColour = "white"
        )
        addStyle(wb, "All_Tested_Pathways", header_style, rows = 1, cols = 1:ncol(summary_df), gridExpand = TRUE)
        addStyle(wb, "Significant_Only", header_style, rows = 1, cols = 1:ncol(summary_df_sig), gridExpand = TRUE)
        
        setColWidths(wb, "All_Tested_Pathways", cols = 1:ncol(summary_df), widths = "auto")
        setColWidths(wb, "Significant_Only", cols = 1:ncol(summary_df_sig), widths = "auto")
        
        saveWorkbook(wb, xlsx_out, overwrite = TRUE)

        cat("  Saved detailed pathway summaries for", coll, "to dataset folder\n")
    }
    cat("\n")
}

cat("====================================================================\n")
cat("GSEA Integration Complete!\n")
cat("====================================================================\n")
