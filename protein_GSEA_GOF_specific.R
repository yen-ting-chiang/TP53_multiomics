################################################################################
# TP53 Protein GSEA - GOF Specific Results Extraction
#
# Reads the Detailed_Pathway_Summary files from the protein_GSEA folder.
# Extracts pathways that meet the following GOF-specific activity criteria:
#   1. Exhibits significant change in GOF_vs_LOF OR Hot_vs_LOF (either Pos or Neg)
#   2. DOES NOT exhibit significant change in LOF_vs_wt (in either Pos or Neg)
#
# This highlights pathways that are uniquely altered by mutant p53 GOF,
# beyond simple loss of function.
#
# Outputs: .csv and .xlsx files saved into protein_GSEA folder.
################################################################################

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(openxlsx)
})

base_path <- "C:/Users/danny/Documents/R_project/TP53_multiomics"
gsea_dir <- file.path(base_path, "protein_GSEA")

datasets <- c("brca_cptac_2020", "coad_cptac_2019", "gbm_cptac_2021",
               "luad_cptac_2020", "paad_cptac_2021", "ucec_cptac_2020")

cat("====================================================================\n")
cat("TP53 Protein GSEA: GOF-Specific Extraction\n")
cat("====================================================================\n\n")

if (!dir.exists(gsea_dir)) {
    stop("GSEA directory not found. Please make sure protein_GSEA exists.")
}

global_gof_list <- list()

for (ds_folder in datasets) {
    ds_out <- file.path(gsea_dir, ds_folder)
    
    if (!dir.exists(ds_out)) {
        next
    }
    
    cat("Extracting specific GOF signals for Dataset:", ds_folder, "\n")

    collections <- list.dirs(ds_out, recursive = FALSE, full.names = FALSE)

    for (coll in collections) {
        coll_dir <- file.path(ds_out, coll)
        if (!dir.exists(coll_dir)) {
            next
        }

        summary_file <- file.path(coll_dir, paste0("Detailed_Pathway_Summary_", ds_folder, "_", coll, ".csv"))
        
        if (!file.exists(summary_file)) {
            next
        }

        # Read the summary dataframe
        df <- fread(summary_file)
        
        # Ensure the required columns exist (sometimes they are NA if no significance found)
        if (!"Pos_Significant_In" %in% names(df)) df$Pos_Significant_In <- ""
        if (!"Neg_Significant_In" %in% names(df)) df$Neg_Significant_In <- ""
        
        # Replace NA with empty string for grepl
        df$Pos_Significant_In[is.na(df$Pos_Significant_In)] <- ""
        df$Neg_Significant_In[is.na(df$Neg_Significant_In)] <- ""
        
        # Condition 1: MUST CONTAIN "GOF_vs_LOF" OR "Hot_vs_LOF" in Pos or Neg
        has_gof_or_hot <- 
            grepl("GOF_vs_LOF", df$Pos_Significant_In) | grepl("GOF_vs_LOF", df$Neg_Significant_In) |
            grepl("Hot_vs_LOF", df$Pos_Significant_In)   | grepl("Hot_vs_LOF",   df$Neg_Significant_In)
            
        # Condition 2: MUST NOT CONTAIN "LOF_vs_wt" in Pos or Neg
        has_lof_wt <-
            grepl("LOF_vs_wt", df$Pos_Significant_In) | grepl("LOF_vs_wt", df$Neg_Significant_In)
            
        # Filter dataframe
        filtered_df <- df %>% filter(has_gof_or_hot & !has_lof_wt)
        
        if (nrow(filtered_df) == 0) {
            cat("  [SKIP]", coll, "- No pathways met the GOF-specific criteria.\n")
            next
        }

        # Output specific file naming
        csv_out <- file.path(coll_dir, paste0("GOF_Specific_Pathway_", ds_folder, "_", coll, ".csv"))
        xlsx_out <- file.path(coll_dir, paste0("GOF_Specific_Pathway_", ds_folder, "_", coll, ".xlsx"))
        
        fwrite(filtered_df, csv_out)
        
        wb <- createWorkbook()
        addWorksheet(wb, "GOF_Specific_Only")
        writeData(wb, "GOF_Specific_Only", filtered_df)
        
        header_style <- createStyle(
            textDecoration = "bold", halign = "center",
            border = "bottom", fgFill = "#C00000", fontColour = "white"
        )
        addStyle(wb, "GOF_Specific_Only", header_style, rows = 1, cols = 1:ncol(filtered_df), gridExpand = TRUE)
        setColWidths(wb, "GOF_Specific_Only", cols = 1:ncol(filtered_df), widths = "auto")
        
        saveWorkbook(wb, xlsx_out, overwrite = TRUE)
        
        cat("  ->", coll, ": Found", nrow(filtered_df), "GOF-specific pathways.\n")
        
        # Add dataset name and store it for global aggregation
        filtered_df$dataset <- ds_folder
        
        if (is.null(global_gof_list[[coll]])) {
            global_gof_list[[coll]] <- list()
        }
        global_gof_list[[coll]][[length(global_gof_list[[coll]]) + 1]] <- filtered_df
    }
    cat("\n")
}

cat("====================================================================\n")
cat("Generating Cross-Dataset GOF-Specific Summaries\n")
cat("====================================================================\n\n")

for (coll in names(global_gof_list)) {
    if (length(global_gof_list[[coll]]) == 0) next
    
    # Bind all GOF specific pathways for this collection across all datasets
    all_gof_df <- bind_rows(global_gof_list[[coll]])
    
    # Calculate pos/neg status per dataset row
    all_gof_df <- all_gof_df %>%
        mutate(
            Is_Pos = grepl("GOF_vs_LOF|Hot_vs_LOF", Pos_Significant_In),
            Is_Neg = grepl("GOF_vs_LOF|Hot_vs_LOF", Neg_Significant_In)
        )
        
    # Summarize across datasets per pathway
    global_summary <- all_gof_df %>%
        group_by(pathway) %>%
        summarize(
            Total_Dataset_Count = n(),
            GOF_Pos_Count = sum(Is_Pos),
            GOF_Neg_Count = sum(Is_Neg),
            Pos_In_Datasets = paste(dataset[Is_Pos], collapse = ", "),
            Neg_In_Datasets = paste(dataset[Is_Neg], collapse = ", ")
        ) %>%
        arrange(desc(Total_Dataset_Count), desc(GOF_Pos_Count), desc(GOF_Neg_Count), pathway)
        
    # Write output
    csv_out <- file.path(gsea_dir, paste0("CrossDataset_GOF_Specific_Summary_", coll, ".csv"))
    xlsx_out <- file.path(gsea_dir, paste0("CrossDataset_GOF_Specific_Summary_", coll, ".xlsx"))
    
    fwrite(global_summary, csv_out)
    
    wb <- createWorkbook()
    addWorksheet(wb, "Global_GOF_Summary")
    writeData(wb, "Global_GOF_Summary", global_summary)
    
    header_style <- createStyle(
        textDecoration = "bold", halign = "center",
        border = "bottom", fgFill = "#7030A0", fontColour = "white"
    )
    addStyle(wb, "Global_GOF_Summary", header_style, rows = 1, cols = 1:ncol(global_summary), gridExpand = TRUE)
    setColWidths(wb, "Global_GOF_Summary", cols = 1:ncol(global_summary), widths = "auto")
    
    saveWorkbook(wb, xlsx_out, overwrite = TRUE)
    
    cat("  Saved Cross-Dataset Summary for:", coll, "\n")
}

cat("\n====================================================================\n")
cat("GOF Specific Extraction Complete!\n")
cat("====================================================================\n")
