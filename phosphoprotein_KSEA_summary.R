################################################################################
# TP53 Phosphoprotein KSEA Results Summarization
#
# Aggregates the results from the KSEA analysis (phosphoprotein_KSEA)
#
# For each dataset across the 9 comparisons, this script tallies for each Kinase.Gene:
#   - Total number of times it is significant (FDR < 0.05)
#   - Total positive significant (FDR < 0.05 & z.score > 0)
#   - Total negative significant (FDR < 0.05 & z.score < 0)
#   - A list of the specific comparisons where it was positive significant.
#   - A list of the specific comparisons where it was negative significant.
#
# Outputs: .csv and .xlsx files saved into phosphoprotein_KSEA folder.
################################################################################

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(openxlsx)
})

base_path <- "C:/Users/danny/Documents/R_project/TP53_multiomics"
ksea_dir <- file.path(base_path, "phosphoprotein_KSEA")

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

cat("====================================================================\n")
cat("TP53 Phosphoprotein KSEA Summarization\n")
cat("====================================================================\n\n")

if (!dir.exists(ksea_dir)) {
    stop("KSEA directory not found. Please run the KSEA script first.")
}

for (i in seq_len(nrow(datasets))) {
    ds_folder <- datasets$folder[i]
    ds_cancer <- datasets$cancer_type[i]

    cat("Processing Summary for Dataset:", ds_folder, "(", ds_cancer, ")\n")

    ds_out <- file.path(ksea_dir, ds_folder)
    if (!dir.exists(ds_out)) {
        cat("  [SKIP] Dataset KSEA directory not found\n")
        next
    }

    # Gather data across all 9 comparisons for this dataset
    all_kinases <- list()

    for (j in seq_len(nrow(comparisons))) {
        comp_name <- comparisons$name[j]
        comp_label <- comparisons$label[j]

        ksea_file <- file.path(ds_out, paste0("KSEA_", comp_name, ".csv"))

        if (file.exists(ksea_file)) {
            ksea_res <- tryCatch(fread(ksea_file), error = function(e) NULL)
            
            if (!is.null(ksea_res) && nrow(ksea_res) > 0 && "Kinase.Gene" %in% colnames(ksea_res)) {
                # Format specific df for binding
                df <- data.frame(
                    Kinase.Gene = ksea_res$Kinase.Gene,
                    FDR = ksea_res$FDR,
                    z.score = ksea_res$z.score,
                    comparison = comp_label,
                    stringsAsFactors = FALSE
                )
                all_kinases[[comp_name]] <- df
            }
        }
    }

    if (length(all_kinases) == 0) {
        cat("  [SKIP] No KSEA results found for this dataset.\n")
        next
    }

    # Bind all rows
    combined_df <- bind_rows(all_kinases)

    # Summarize per Kinase.Gene
    summary_df <- combined_df %>%
        group_by(Kinase.Gene) %>%
        summarize(
            Total_Tested = n(),
            Sig_Count = sum(!is.na(FDR) & FDR < 0.05),
            Pos_Sig_Count = sum(!is.na(FDR) & FDR < 0.05 & !is.na(z.score) & z.score > 0),
            Neg_Sig_Count = sum(!is.na(FDR) & FDR < 0.05 & !is.na(z.score) & z.score < 0),
            Pos_Significant_In = paste(comparison[!is.na(FDR) & FDR < 0.05 & !is.na(z.score) & z.score > 0], collapse = ", "),
            Neg_Significant_In = paste(comparison[!is.na(FDR) & FDR < 0.05 & !is.na(z.score) & z.score < 0], collapse = ", ")
        ) %>%
        arrange(desc(Sig_Count), Kinase.Gene)

    # Filter to only keep kinases that are strictly significant at least once
    summary_df_sig <- summary_df %>% filter(Sig_Count > 0)

    # Write per-dataset outputs
    csv_out <- file.path(ds_out, paste0("Detailed_Kinase_Summary_", ds_folder, ".csv"))
    fwrite(summary_df, csv_out)
    
    xlsx_out <- file.path(ds_out, paste0("Detailed_Kinase_Summary_", ds_folder, ".xlsx"))
    wb <- createWorkbook()
    
    addWorksheet(wb, "All_Tested_Kinases")
    writeData(wb, "All_Tested_Kinases", summary_df)
    
    addWorksheet(wb, "Significant_Only")
    writeData(wb, "Significant_Only", summary_df_sig)
    
    # Styles for Excel
    header_style <- createStyle(
        textDecoration = "bold", halign = "center",
        border = "bottom", fgFill = "#4472C4", fontColour = "white"
    )
    addStyle(wb, "All_Tested_Kinases", header_style, rows = 1, cols = 1:ncol(summary_df), gridExpand = TRUE)
    addStyle(wb, "Significant_Only", header_style, rows = 1, cols = 1:ncol(summary_df_sig), gridExpand = TRUE)
    
    setColWidths(wb, "All_Tested_Kinases", cols = 1:ncol(summary_df), widths = "auto")
    setColWidths(wb, "Significant_Only", cols = 1:ncol(summary_df_sig), widths = "auto")
    
    saveWorkbook(wb, xlsx_out, overwrite = TRUE)

    cat("  Saved detailed summaries to dataset folder\n\n")
}

cat("====================================================================\n")
cat("KSEA Integration Complete!\n")
cat("====================================================================\n")
