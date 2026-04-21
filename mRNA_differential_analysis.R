################################################################################
# TP53 mRNA-Level Differential Expression Analysis (limma)
#
# Comparisons:
#   1. TP53mt vs TP53wt
#   2. MUT_GOF vs MUT_LOF
#   3. Hotspots vs MUT_LOF
#   4. MUT_GOF vs TP53wt
#   5. MUT_LOF vs TP53wt
#   6. Hotspots vs TP53wt
#   7. DN vs TP53wt
#   8. Non-DN vs TP53wt
#   9. DN vs non-DN
#
# Covariates: purity, sex, age (NO batch, following mRNA_GSEA_analysis pipeline)
# GBM dataset: log2-transform before limma; others: direct limma
#
# Methodology adapted from mRNA_GSEA_analysis/ pipeline
################################################################################

# ==============================================================================
# Section 1: Setup
# ==============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(janitor)
    library(glue)
    library(limma)
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

# Output directory
output_dir <- file.path(base_path, "mRNA_differential_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Minimum fraction of non-NA values per gene
min_frac_complete <- 0.75
# Minimum samples per group
min_per_group <- 5

cat("====================================================================\n")
cat("TP53 mRNA Differential Expression Analysis\n")
cat("====================================================================\n\n")

# ==============================================================================
# Section 2: Utility Functions (adapted from mRNA_GSEA_analysis)
# ==============================================================================

#' Read case list from cBioPortal format file
read_case_list <- function(path_file) {
    if (!file.exists(path_file)) return(character(0))
    x <- readLines(path_file, warn = FALSE, encoding = "UTF-8")
    line <- x[grepl("^case_list_ids:", x)]
    if (!length(line)) return(character(0))
    ids <- sub("^case_list_ids:\\s*", "", line[1])
    ids <- unlist(strsplit(ids, "[,\\s]+"))
    unique(ids[nchar(ids) > 0])
}

#' Load mRNA expression matrix from dataset directory
#' Tries multiple filenames (fpkm, rsem, rpkm) in order
load_mrna_matrix <- function(ds_dir) {
    candidates <- c(
        "data_mrna_seq_fpkm.txt",
        "data_mrna_seq_v2_rsem.txt",
        "data_mrna_seq_rpkm.txt",
        "data_mrna_seq_rsem.txt"
    )
    hit_flags <- file.exists(file.path(ds_dir, candidates))
    if (!any(hit_flags)) {
        stop(paste("RNA file not found in:", ds_dir,
                   "\nTried:", paste(candidates, collapse = ", ")))
    }
    fp <- file.path(ds_dir, candidates[which(hit_flags)[1]])
    cat("  Reading mRNA matrix:", basename(fp), "\n")

    dat <- suppressMessages(readr::read_tsv(fp, guess_max = 200000, show_col_types = FALSE))

    # Identify gene column
    gene_cols <- c("Hugo_Symbol", "hugo_symbol", "Gene", "Gene_Symbol",
                   "HugoSymbol", "GENE_SYMBOL", "gene", "gene_symbol")
    gcol <- intersect(gene_cols, names(dat))
    if (!length(gcol)) gcol <- names(dat)[1]
    dat <- dplyr::rename(dat, Gene = !!gcol[1])
    dat$Gene <- sub("\\|.*$", "", dat$Gene)

    # Identify sample columns
    not_sample <- c("Gene", "Entrez_Gene_Id", "Entrez_Gene_Id.",
                    "ENTREZ_GENE_ID", "Description", "Gene_Name",
                    "GeneName", "Gene_Symbol")
    sample_cols_all <- setdiff(names(dat), not_sample)

    # RNA-specific case list
    bn <- basename(fp)
    case_candidates <- c(
        file.path(ds_dir, "case_lists", sub("^data_", "cases_", bn)),
        file.path(ds_dir, "case_lists", "cases_mrna_seq_fpkm.txt"),
        file.path(ds_dir, "case_lists", "cases_mrna_seq_v2_rsem.txt"),
        file.path(ds_dir, "case_lists", "cases_mrna_seq_rpkm.txt"),
        file.path(ds_dir, "case_lists", "cases_mrna_seq_rsem.txt")
    )
    case_hit <- case_candidates[file.exists(case_candidates)]
    keep_ids <- if (length(case_hit)) read_case_list(case_hit[1]) else character(0)

    if (length(keep_ids)) {
        inter <- intersect(sample_cols_all, keep_ids)
        sample_cols <- if (length(inter) >= 10) inter else sample_cols_all
    } else {
        sample_cols <- sample_cols_all
    }

    # Create matrix
    m <- dat %>%
        dplyr::select(Gene, dplyr::all_of(sample_cols)) %>%
        janitor::remove_empty("cols")

    rn <- as.character(m$Gene)
    ok <- !is.na(rn) & nzchar(trimws(rn))
    m <- m[ok, , drop = FALSE]
    rn <- rn[ok]
    m <- as.matrix(m[, -1, drop = FALSE])
    storage.mode(m) <- "double"
    rownames(m) <- rn

    # Handle duplicates
    if (anyDuplicated(rownames(m))) {
        cat("  Averaging duplicate genes\n")
        m <- rowsum(m, group = rownames(m), reorder = FALSE) /
            as.vector(table(rownames(m)))
    }

    cat("  Matrix dimensions:", nrow(m), "genes x", ncol(m), "samples\n")
    m
}

#' Impute and filter (gene-median imputation, following mRNA pipeline)
impute_and_filter <- function(mat, min_frac = 0.75) {
    keep <- rowMeans(is.finite(mat)) >= min_frac
    m <- mat[keep, , drop = FALSE]
    if (anyNA(m)) {
        meds <- apply(m, 1, function(v) median(v[is.finite(v)], na.rm = TRUE))
        for (i in seq_len(nrow(m))) {
            vi <- m[i, ]
            vi[!is.finite(vi)] <- meds[i]
            m[i, ] <- vi
        }
    }
    m
}

#' Get tumor purity covariate (from clinical patient file)
get_purity_covariate <- function(ds_dir, sample_ids) {
    fp <- file.path(ds_dir, "data_clinical_patient.txt")
    if (!file.exists(fp)) return(NULL)

    clin <- tryCatch(
        readr::read_tsv(fp, comment = "#", show_col_types = FALSE),
        error = function(e) NULL
    )
    if (is.null(clin)) return(NULL)

    purity_cols <- c("TUMOR_PURITY", "Tumor_Purity", "tumor_purity",
                     "PURITY", "Purity", "purity")
    hit <- intersect(purity_cols, colnames(clin))
    if (!length(hit)) return(NULL)

    patient_id_cols <- c("PATIENT_ID", "Patient_ID", "patient_id")
    pid_col <- intersect(patient_id_cols, colnames(clin))
    if (!length(pid_col)) return(NULL)

    clin_df <- data.frame(
        patient_id = as.character(clin[[pid_col[1]]]),
        purity = suppressWarnings(as.numeric(clin[[hit[1]]])),
        stringsAsFactors = FALSE
    )

    purity_vec <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
    for (i in seq_along(sample_ids)) {
        sid <- sample_ids[i]
        matches <- clin_df$patient_id[clin_df$patient_id == sid |
                                          sapply(clin_df$patient_id, function(pid) grepl(pid, sid, fixed = TRUE))]
        if (length(matches) > 0) {
            purity_vec[sid] <- clin_df$purity[clin_df$patient_id == matches[1]]
        }
    }
    purity_vec
}

#' Get sex and age covariates (from clinical sample file)
get_sex_age_covariates <- function(ds_dir, sample_ids) {
    fp <- file.path(ds_dir, "data_clinical_sample.txt")
    empty_df <- data.frame(
        sex = rep(NA_character_, length(sample_ids)),
        age = rep(NA_real_, length(sample_ids)),
        row.names = sample_ids, stringsAsFactors = FALSE
    )
    if (!file.exists(fp)) return(empty_df)

    clin <- tryCatch(
        readr::read_tsv(fp, comment = "#", show_col_types = FALSE),
        error = function(e) NULL
    )
    if (is.null(clin)) return(empty_df)

    id_cols <- intersect(c("SAMPLE_ID", "sample_id", "Sample_ID"), colnames(clin))
    if (!length(id_cols)) return(empty_df)
    clin$SAMPLE_ID <- as.character(clin[[id_cols[1]]])

    sex_cols <- c("SEX", "Sex", "sex", "GENDER", "Gender", "gender")
    sex_hit <- intersect(sex_cols, colnames(clin))
    sex_values <- if (length(sex_hit)) as.character(clin[[sex_hit[1]]]) else rep(NA_character_, nrow(clin))

    age_cols <- c("AGE", "Age", "age", "AGE_AT_DIAGNOSIS", "Age_at_Diagnosis")
    age_hit <- intersect(age_cols, colnames(clin))
    age_values <- if (length(age_hit)) suppressWarnings(as.numeric(clin[[age_hit[1]]])) else rep(NA_real_, nrow(clin))

    clin_aligned <- data.frame(
        SAMPLE_ID = clin$SAMPLE_ID, sex = sex_values, age = age_values,
        stringsAsFactors = FALSE
    )
    clin_aligned <- clin_aligned[match(sample_ids, clin_aligned$SAMPLE_ID), ]
    rownames(clin_aligned) <- sample_ids
    clin_aligned[, c("sex", "age"), drop = FALSE]
}

#' Safely coerce covariates
coerce_covariates_safely <- function(df) {
    df <- as.data.frame(df, check.names = FALSE)
    keep <- rep(TRUE, ncol(df))
    names(keep) <- colnames(df)
    for (cn in colnames(df)) {
        v <- df[[cn]]
        if (is.factor(v) || is.character(v) || is.logical(v)) {
            v <- factor(v)
            lv <- levels(droplevels(v[!is.na(v)]))
            if (length(lv) <= 1) {
                keep[cn] <- FALSE
                cat("    [covars] drop single-level factor:", cn, "\n")
            } else {
                df[[cn]] <- v
            }
        } else {
            df[[cn]] <- suppressWarnings(as.numeric(v))
        }
    }
    df[, keep, drop = FALSE]
}

# ==============================================================================
# Section 3: Limma DEG Function for Group Comparison
# ==============================================================================

#' Run limma DEG analysis for a binary group comparison
#' Covariates: purity, sex, age (NO batch, following mRNA_GSEA_analysis)
run_limma_group_comparison <- function(mat, group_vec,
                                       purity_vec = NULL,
                                       sa_df = NULL) {
    # Align samples
    common <- intersect(colnames(mat), names(group_vec))
    group_vec <- group_vec[common]
    group_vec <- droplevels(group_vec)

    if (nlevels(group_vec) < 2) {
        cat("    [SKIP] Less than 2 group levels\n")
        return(NULL)
    }

    tab <- table(group_vec)
    if (any(tab < min_per_group)) {
        cat("    [SKIP] Group too small:", paste(names(tab), "=", tab, collapse = ", "), "\n")
        return(NULL)
    }

    M <- mat[, common, drop = FALSE]
    so <- common

    # Build covariate data frame (NO batch, following mRNA pipeline)
    DF <- data.frame(group = group_vec, row.names = so)

    # Add purity
    if (!is.null(purity_vec)) {
        p <- suppressWarnings(as.numeric(purity_vec[so]))
        if (sum(is.finite(p)) >= length(so) * 0.6) {
            DF$purity <- p
        }
    }

    # Add sex and age
    if (!is.null(sa_df)) {
        if ("sex" %in% colnames(sa_df)) {
            s <- sa_df[so, "sex"]
            if (sum(!is.na(s)) >= length(so) * 0.8) {
                DF$sex <- factor(s)
            }
        }
        if ("age" %in% colnames(sa_df)) {
            a <- suppressWarnings(as.numeric(sa_df[so, "age"]))
            if (sum(is.finite(a)) >= length(so) * 0.8) {
                DF$age <- a
            }
        }
    }

    # Remove all-NA columns
    all_na <- vapply(DF, function(z) all(is.na(z)), logical(1))
    if (any(all_na)) DF <- DF[, !all_na, drop = FALSE]

    # Clean covariates
    DF <- coerce_covariates_safely(DF)

    if (!("group" %in% colnames(DF))) {
        cat("    [ERROR] Group variable dropped\n")
        return(NULL)
    }

    ok <- !is.na(DF$group)
    DF <- DF[ok, , drop = FALSE]
    M <- M[, rownames(DF), drop = FALSE]

    # Build design matrix
    des <- model.matrix(~ 0 + ., data = DF)

    lvls <- levels(DF$group)
    group_cols <- paste0("group", lvls)
    if (!all(group_cols %in% colnames(des))) {
        cat("    [ERROR] Group columns not found in design matrix\n")
        return(NULL)
    }

    contrast_str <- paste0(group_cols[2], " - ", group_cols[1])
    contrast_mat <- limma::makeContrasts(contrasts = contrast_str, levels = des)

    # Filter genes with insufficient df
    rnk <- qr(des)$rank
    need <- rnk + 1L
    nobs <- rowSums(is.finite(M))
    keep_rows <- nobs >= need
    if (sum(keep_rows) == 0) {
        cat("    [SKIP] No genes pass df filter\n")
        return(NULL)
    }
    M <- M[keep_rows, , drop = FALSE]

    # Run limma
    fit <- limma::lmFit(M, design = des)
    fit2 <- limma::contrasts.fit(fit, contrast_mat)
    eb <- limma::eBayes(fit2, trend = TRUE)

    tbl <- limma::topTable(eb, number = Inf, sort.by = "P")
    tbl$gene <- rownames(tbl)
    tbl <- tbl[, c("gene", "logFC", "t", "P.Value", "adj.P.Val", "B")]

    cat("    Results:", nrow(tbl), "genes tested,",
        sum(tbl$adj.P.Val < 0.05, na.rm = TRUE), "significant (padj < 0.05)\n")
    tbl
}

# ==============================================================================
# Section 4: Helper to run one comparison and save results
# ==============================================================================

run_and_save_comparison <- function(comparison_name, file_prefix, sheet_name,
                                    group_vec, mat, purity_vec, sa_df,
                                    ds_out, ds_summary) {
    cat("\n  ---", comparison_name, "---\n")

    deg <- run_limma_group_comparison(mat, group_vec, purity_vec, sa_df)

    prefix_clean <- gsub("[^A-Za-z0-9_]", "_", file_prefix)

    if (!is.null(deg)) {
        fwrite(deg, file.path(ds_out, paste0("DEG_", prefix_clean, ".csv")))
        wb <- createWorkbook()
        addWorksheet(wb, sheet_name)
        writeData(wb, 1, deg)
        saveWorkbook(wb, file.path(ds_out, paste0("DEG_", prefix_clean, ".xlsx")), overwrite = TRUE)

        ds_summary[[paste0(prefix_clean, "_total")]] <- nrow(deg)
        ds_summary[[paste0(prefix_clean, "_sig")]] <- sum(deg$adj.P.Val < 0.05, na.rm = TRUE)
        ds_summary[[paste0(prefix_clean, "_up")]] <- sum(deg$adj.P.Val < 0.05 & deg$logFC > 0, na.rm = TRUE)
        ds_summary[[paste0(prefix_clean, "_down")]] <- sum(deg$adj.P.Val < 0.05 & deg$logFC < 0, na.rm = TRUE)
    } else {
        ds_summary[[paste0(prefix_clean, "_total")]] <- NA
        ds_summary[[paste0(prefix_clean, "_sig")]] <- NA
        ds_summary[[paste0(prefix_clean, "_up")]] <- NA
        ds_summary[[paste0(prefix_clean, "_down")]] <- NA
    }
    ds_summary
}

# ==============================================================================
# Section 5: Main Processing Loop
# ==============================================================================

# Load TP53 classification results
tp53_file <- file.path(base_path, "TP53_mutation_classification",
                       "all_CPTAC_TP53_classification.csv")
tp53_all <- read.csv(tp53_file, stringsAsFactors = FALSE)
cat("Loaded TP53 classification:", nrow(tp53_all), "samples\n\n")

all_deg_summary <- list()

for (i in seq_len(nrow(datasets))) {
    ds_folder <- datasets$folder[i]
    ds_cancer <- datasets$cancer_type[i]
    ds_dir <- file.path(base_path, ds_folder)

    cat("==================================================================\n")
    cat("Processing:", ds_folder, "(", ds_cancer, ")\n")
    cat("==================================================================\n")

    # --- Load mRNA matrix ---
    mat0 <- load_mrna_matrix(ds_dir)

    # Log-transform: ONLY for GBM (per user instruction)
    if (ds_folder == "gbm_cptac_2021") {
        cat("  Log2-transforming GBM data\n")
        mat0 <- log2(mat0 + 1)
    }

    # Impute and filter
    mat <- impute_and_filter(mat0, min_frac = min_frac_complete)
    cat("  After filtering:", nrow(mat), "genes x", ncol(mat), "samples\n")

    # --- Get covariates (NO batch per mRNA pipeline) ---
    purity_vec <- get_purity_covariate(ds_dir, colnames(mat))
    sa_df <- get_sex_age_covariates(ds_dir, colnames(mat))
    cat("  Covariates loaded (purity, sex, age; no batch per mRNA pipeline)\n")

    # --- Get TP53 classification for this dataset ---
    tp53_ds <- tp53_all %>% filter(cancer_type == ds_cancer)
    tp53_mut <- tp53_ds %>% filter(mt == 1)

    # Pre-define sample groups
    wt_samples <- tp53_ds$sample_id[tp53_ds$wt == 1]
    gof_samples <- tp53_mut$sample_id[tp53_mut$GOF == 1]
    lof_samples <- tp53_mut$sample_id[tp53_mut$LOF == 1]
    hotspot_samples <- tp53_mut$sample_id[tp53_mut$hotspot == 1]
    dn_samples <- tp53_mut$sample_id[tp53_mut$DN == 1]
    nondn_samples <- tp53_mut$sample_id[tp53_mut$non_DN == 1]

    # Create output subdirectory
    ds_out <- file.path(output_dir, ds_folder)
    dir.create(ds_out, recursive = TRUE, showWarnings = FALSE)

    ds_summary <- list(dataset = ds_folder, cancer_type = ds_cancer)

    # ========================================================================
    # Comparison 1: TP53mt vs TP53wt
    # ========================================================================
    g1 <- setNames(
        ifelse(tp53_ds$mt == 1, "TP53mt", "TP53wt"), tp53_ds$sample_id
    )
    g1 <- factor(g1[!is.na(g1)], levels = c("TP53wt", "TP53mt"))
    cat("    Group sizes: TP53wt=", sum(g1 == "TP53wt"), " TP53mt=", sum(g1 == "TP53mt"), "\n")
    ds_summary <- run_and_save_comparison(
        "Comparison 1: TP53mt vs TP53wt", "TP53mt_vs_TP53wt", "mt_vs_wt",
        g1, mat, purity_vec, sa_df, ds_out, ds_summary
    )

    # ========================================================================
    # Comparison 2: MUT_GOF vs MUT_LOF
    # ========================================================================
    g2 <- setNames(
        c(rep("MUT_GOF", length(gof_samples)), rep("MUT_LOF", length(lof_samples))),
        c(gof_samples, lof_samples)
    )
    g2 <- factor(g2, levels = c("MUT_LOF", "MUT_GOF"))
    cat("    GOF:", length(gof_samples), " LOF:", length(lof_samples), "\n")
    ds_summary <- run_and_save_comparison(
        "Comparison 2: MUT_GOF vs MUT_LOF", "MUT_GOF_vs_MUT_LOF", "GOF_vs_LOF",
        g2, mat, purity_vec, sa_df, ds_out, ds_summary
    )

    # ========================================================================
    # Comparison 3: Hotspots vs MUT_LOF
    # ========================================================================
    g3_labels <- setNames(
        c(rep("Hotspot", length(hotspot_samples)), rep("MUT_LOF", length(lof_samples))),
        c(hotspot_samples, lof_samples)
    )
    g3_labels <- g3_labels[!duplicated(names(g3_labels))]
    g3 <- factor(g3_labels, levels = c("MUT_LOF", "Hotspot"))
    cat("    Hotspot:", length(hotspot_samples), " LOF:", length(lof_samples), "\n")
    ds_summary <- run_and_save_comparison(
        "Comparison 3: Hotspot vs MUT_LOF", "Hotspot_vs_MUT_LOF", "Hot_vs_LOF",
        g3, mat, purity_vec, sa_df, ds_out, ds_summary
    )

    # ========================================================================
    # Comparison 4: MUT_GOF vs TP53wt
    # ========================================================================
    g4 <- setNames(
        c(rep("MUT_GOF", length(gof_samples)), rep("TP53wt", length(wt_samples))),
        c(gof_samples, wt_samples)
    )
    g4 <- factor(g4, levels = c("TP53wt", "MUT_GOF"))
    cat("    GOF:", length(gof_samples), " TP53wt:", length(wt_samples), "\n")
    ds_summary <- run_and_save_comparison(
        "Comparison 4: MUT_GOF vs TP53wt", "MUT_GOF_vs_TP53wt", "GOF_vs_wt",
        g4, mat, purity_vec, sa_df, ds_out, ds_summary
    )

    # ========================================================================
    # Comparison 5: MUT_LOF vs TP53wt
    # ========================================================================
    g5 <- setNames(
        c(rep("MUT_LOF", length(lof_samples)), rep("TP53wt", length(wt_samples))),
        c(lof_samples, wt_samples)
    )
    g5 <- factor(g5, levels = c("TP53wt", "MUT_LOF"))
    cat("    LOF:", length(lof_samples), " TP53wt:", length(wt_samples), "\n")
    ds_summary <- run_and_save_comparison(
        "Comparison 5: MUT_LOF vs TP53wt", "MUT_LOF_vs_TP53wt", "LOF_vs_wt",
        g5, mat, purity_vec, sa_df, ds_out, ds_summary
    )

    # ========================================================================
    # Comparison 6: Hotspots vs TP53wt
    # ========================================================================
    g6 <- setNames(
        c(rep("Hotspot", length(hotspot_samples)), rep("TP53wt", length(wt_samples))),
        c(hotspot_samples, wt_samples)
    )
    g6 <- factor(g6, levels = c("TP53wt", "Hotspot"))
    cat("    Hotspot:", length(hotspot_samples), " TP53wt:", length(wt_samples), "\n")
    ds_summary <- run_and_save_comparison(
        "Comparison 6: Hotspot vs TP53wt", "Hotspot_vs_TP53wt", "Hot_vs_wt",
        g6, mat, purity_vec, sa_df, ds_out, ds_summary
    )

    # ========================================================================
    # Comparison 7: DN vs TP53wt
    # ========================================================================
    g7 <- setNames(
        c(rep("DN", length(dn_samples)), rep("TP53wt", length(wt_samples))),
        c(dn_samples, wt_samples)
    )
    g7 <- factor(g7, levels = c("TP53wt", "DN"))
    cat("    DN:", length(dn_samples), " TP53wt:", length(wt_samples), "\n")
    ds_summary <- run_and_save_comparison(
        "Comparison 7: DN vs TP53wt", "DN_vs_TP53wt", "DN_vs_wt",
        g7, mat, purity_vec, sa_df, ds_out, ds_summary
    )

    # ========================================================================
    # Comparison 8: Non-DN vs TP53wt
    # ========================================================================
    g8 <- setNames(
        c(rep("NonDN", length(nondn_samples)), rep("TP53wt", length(wt_samples))),
        c(nondn_samples, wt_samples)
    )
    g8 <- factor(g8, levels = c("TP53wt", "NonDN"))
    cat("    NonDN:", length(nondn_samples), " TP53wt:", length(wt_samples), "\n")
    ds_summary <- run_and_save_comparison(
        "Comparison 8: Non-DN vs TP53wt", "NonDN_vs_TP53wt", "NonDN_vs_wt",
        g8, mat, purity_vec, sa_df, ds_out, ds_summary
    )

    # ========================================================================
    # Comparison 9: DN vs non-DN
    # ========================================================================
    g9_labels <- setNames(
        c(rep("DN", length(dn_samples)), rep("NonDN", length(nondn_samples))),
        c(dn_samples, nondn_samples)
    )
    g9_labels <- g9_labels[!duplicated(names(g9_labels))]
    g9 <- factor(g9_labels, levels = c("NonDN", "DN"))
    cat("    DN:", length(dn_samples), " NonDN:", length(nondn_samples), "\n")
    ds_summary <- run_and_save_comparison(
        "Comparison 9: DN vs non-DN", "DN_vs_NonDN", "DN_vs_NonDN",
        g9, mat, purity_vec, sa_df, ds_out, ds_summary
    )

    all_deg_summary[[ds_folder]] <- as.data.frame(ds_summary, stringsAsFactors = FALSE)
    cat("\n")
}

# ==============================================================================
# Section 6: Summary Statistics
# ==============================================================================

cat("====================================================================\n")
cat("Generating Summary Statistics\n")
cat("====================================================================\n\n")

summary_df <- bind_rows(all_deg_summary)

# Save summary CSV
fwrite(summary_df, file.path(output_dir, "DEG_summary_statistics.csv"))

# Save summary XLSX with formatting
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
saveWorkbook(wb_sum, file.path(output_dir, "DEG_summary_statistics.xlsx"), overwrite = TRUE)

# Print summary
cat("\n=== mRNA Differential Gene Summary (padj < 0.05) ===\n\n")
print(as.data.frame(summary_df), row.names = FALSE)

cat("\n====================================================================\n")
cat("Pipeline Complete!\n")
cat("All results saved to:", output_dir, "\n")
cat("====================================================================\n")
