################################################################################
# TP53 Protein-Level Differential Expression Analysis (limma)
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
# Covariates: TMT plex (batch), sex, age, tumor purity
#
# Reference: Shahbandi et al. (2023) Cell Death Discovery
# Methodology adapted from protein_level_DEG_analysis/ pipeline
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
    library(imputeLCMD)
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
output_dir <- file.path(base_path, "protein_differential_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Minimum fraction of non-NA values per gene
min_frac_complete <- 0.75
# Minimum samples per group
min_per_group <- 5

cat("====================================================================\n")
cat("TP53 Protein Differential Expression Analysis\n")
cat("====================================================================\n\n")

# ==============================================================================
# Section 2: Utility Functions (adapted from protein_level_DEG_analysis)
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

#' Load protein quantification matrix
load_protein_matrix <- function(ds_dir) {
    fp <- file.path(ds_dir, "data_protein_quantification.txt")
    if (!file.exists(fp)) stop(paste("File not found:", fp))

    cat("  Reading protein matrix:", basename(fp), "\n")
    dat <- suppressMessages(readr::read_tsv(fp, guess_max = 200000, show_col_types = FALSE))

    # Identify gene column
    gene_cols <- c("Hugo_Symbol", "hugo_symbol", "Gene", "Gene_Symbol",
                   "HugoSymbol", "GENE_SYMBOL", "gene", "gene_symbol",
                   "Composite.Element.REF")
    gcol <- intersect(gene_cols, names(dat))
    if (!length(gcol)) gcol <- names(dat)[1]

    dat <- dplyr::rename(dat, Gene = !!gcol[1])
    dat$Gene <- sub("\\|.*$", "", dat$Gene)

    # Identify sample columns
    not_sample <- c("Gene", "Entrez_Gene_Id", "Entrez_Gene_Id.",
                    "ENTREZ_GENE_ID", "Description", "Gene_Name",
                    "GeneName", "Gene_Symbol")
    sample_cols_all <- setdiff(names(dat), not_sample)

    # Filter by case list if available
    case_file <- file.path(ds_dir, "case_lists", "cases_protein_quantification.txt")
    keep_ids <- read_case_list(case_file)

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

    rn <- m$Gene
    m <- as.matrix(m[, -1, drop = FALSE])
    storage.mode(m) <- "double"
    rownames(m) <- rn

    # Handle duplicate genes by averaging
    if (anyDuplicated(rownames(m))) {
        cat("  Averaging duplicate genes\n")
        m <- rowsum(m, group = rownames(m), reorder = FALSE) /
            as.vector(table(rownames(m)))
    }

    cat("  Matrix dimensions:", nrow(m), "genes x", ncol(m), "samples\n")
    m
}

#' Impute missing values and filter low-coverage genes
impute_and_filter <- function(mat, min_frac = 0.75) {
    keep <- rowMeans(!is.na(mat)) >= min_frac
    m <- mat[keep, , drop = FALSE]
    if (any(is.na(m))) {
        m <- imputeLCMD::impute.MinProb(m, q = 0.01)
    }
    m
}

#' Sanitize batch factor levels
sanitize_batch_levels <- function(x, pipe_policy = "NA", min_per_level = 2) {
    x0 <- as.character(x)
    has_pipe <- grepl("\\|", x0)
    if (any(has_pipe, na.rm = TRUE)) {
        x0[has_pipe] <- if (identical(pipe_policy, "NA")) NA_character_ else "b_small"
    }
    fac <- factor(make.names(x0))
    fac <- droplevels(fac)
    if (!is.null(min_per_level) && min_per_level > 1) {
        tab <- table(fac, useNA = "no")
        small <- names(tab)[tab < min_per_level]
        if (length(small)) {
            fac_chr <- as.character(fac)
            fac_chr[fac_chr %in% small] <- "b_small"
            fac <- droplevels(factor(fac_chr))
        }
    }
    fac
}

#' Get batch factor from clinical metadata
get_batch_factor <- function(ds_dir, sample_ids) {
    meta_fp <- file.path(ds_dir, "data_clinical_sample.txt")
    if (!file.exists(meta_fp)) return(NULL)

    meta <- suppressMessages(
        readr::read_tsv(meta_fp, show_col_types = FALSE, comment = "#")
    ) |> as.data.frame()

    id_cols <- intersect(
        c("SAMPLE_ID", "sample_id", "Sample_ID", "Sample", "sample"),
        names(meta)
    )
    if (!length(id_cols)) return(NULL)

    meta <- dplyr::rename(meta, SAMPLE_ID = !!id_cols[1])
    meta$SAMPLE_ID <- as.character(meta$SAMPLE_ID)
    meta <- meta[match(sample_ids, meta$SAMPLE_ID), , drop = FALSE]
    rownames(meta) <- sample_ids

    # Try TMT_PLEX or EXPERIMENT
    cand <- c("TMT_PLEX", "EXPERIMENT")
    hit <- intersect(cand, colnames(meta))

    for (cn in hit) {
        fac <- sanitize_batch_levels(meta[[cn]])
        if (nlevels(fac) >= 2 && sum(!is.na(fac)) >= 3) {
            names(fac) <- sample_ids
            return(list(name = cn, fac = fac))
        }
    }
    NULL
}

#' Get tumor purity covariate
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

#' Get sex and age covariates
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

#' Safely coerce covariates (remove single-level factors, etc.)
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
#'
#' @param mat Expression matrix (genes x samples), already imputed
#' @param group_vec Named factor with two levels (test vs reference)
#' @param batch_fac Batch factor (optional)
#' @param purity_vec Purity numeric vector (optional)
#' @param sa_df Sex/age data frame (optional)
#' @return Data frame with DEG results, or NULL if insufficient data
run_limma_group_comparison <- function(mat, group_vec,
                                       batch_fac = NULL,
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
    so <- common  # sample order

    # Build covariate data frame
    DF <- data.frame(group = group_vec, row.names = so)

    # Add batch
    if (!is.null(batch_fac)) {
        b <- batch_fac[so]
        if (sum(!is.na(b)) >= 3 && nlevels(droplevels(factor(b[!is.na(b)]))) >= 2) {
            DF$batch <- factor(b)
        }
    }

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

    # Ensure group is still present
    if (!("group" %in% colnames(DF))) {
        cat("    [ERROR] Group variable dropped during covariate cleaning\n")
        return(NULL)
    }

    # Keep only complete cases for group
    ok <- !is.na(DF$group)
    DF <- DF[ok, , drop = FALSE]
    M <- M[, rownames(DF), drop = FALSE]

    # Build design matrix
    des <- model.matrix(~ 0 + ., data = DF)

    # Make contrast: test level - reference level (second level vs first)
    lvls <- levels(DF$group)
    # Column names in design matrix
    group_cols <- paste0("group", lvls)

    # Check both group columns exist
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
# Section 4: Main Processing Loop
# ==============================================================================

# Load TP53 classification results
tp53_file <- file.path(base_path, "TP53_mutation_classification",
                       "all_CPTAC_TP53_classification.csv")
tp53_all <- read.csv(tp53_file, stringsAsFactors = FALSE)
cat("Loaded TP53 classification:", nrow(tp53_all), "samples\n\n")

# Storage for summary
all_deg_summary <- list()

for (i in seq_len(nrow(datasets))) {
    ds_folder <- datasets$folder[i]
    ds_cancer <- datasets$cancer_type[i]
    ds_dir <- file.path(base_path, ds_folder)

    cat("==================================================================\n")
    cat("Processing:", ds_folder, "(", ds_cancer, ")\n")
    cat("==================================================================\n")

    # --- Load protein matrix ---
    mat0 <- load_protein_matrix(ds_dir)

    # Log-transform if needed (values > 100 suggest raw counts/intensities)
    mx <- suppressWarnings(max(mat0, na.rm = TRUE))
    if (is.finite(mx) && mx > 100) {
        cat("  Log2-transforming data\n")
        mat0 <- log2(mat0 + 1)
    }

    # Impute and filter
    mat <- impute_and_filter(mat0, min_frac = min_frac_complete)
    cat("  After filtering:", nrow(mat), "genes x", ncol(mat), "samples\n")

    # --- Get covariates ---
    bi <- get_batch_factor(ds_dir, colnames(mat))
    batch_fac <- if (!is.null(bi)) droplevels(bi$fac[colnames(mat)]) else NULL
    if (!is.null(bi)) cat("  Batch variable:", bi$name, "with", nlevels(batch_fac), "levels\n")

    purity_vec <- get_purity_covariate(ds_dir, colnames(mat))
    sa_df <- get_sex_age_covariates(ds_dir, colnames(mat))
    cat("  Covariates loaded (purity, sex, age)\n")

    # --- Get TP53 classification for this dataset ---
    tp53_ds <- tp53_all %>% filter(cancer_type == ds_cancer)

    # Create output subdirectory
    ds_out <- file.path(output_dir, ds_folder)
    dir.create(ds_out, recursive = TRUE, showWarnings = FALSE)

    # Track summary for this dataset
    ds_summary <- data.frame(
        dataset = ds_folder,
        cancer_type = ds_cancer,
        stringsAsFactors = FALSE
    )

    # ========================================================================
    # Comparison 1: TP53mt vs TP53wt
    # ========================================================================
    cat("\n  --- Comparison 1: TP53mt vs TP53wt ---\n")

    group1 <- setNames(
        ifelse(tp53_ds$mt == 1, "TP53mt", "TP53wt"),
        tp53_ds$sample_id
    )
    group1 <- factor(group1[!is.na(group1)], levels = c("TP53wt", "TP53mt"))

    cat("    Group sizes: TP53wt=", sum(group1 == "TP53wt", na.rm = TRUE),
        " TP53mt=", sum(group1 == "TP53mt", na.rm = TRUE), "\n")

    deg1 <- run_limma_group_comparison(mat, group1, batch_fac, purity_vec, sa_df)

    if (!is.null(deg1)) {
        fwrite(deg1, file.path(ds_out, "DEG_TP53mt_vs_TP53wt.csv"))

        wb <- createWorkbook()
        addWorksheet(wb, "TP53mt_vs_TP53wt")
        writeData(wb, 1, deg1)
        saveWorkbook(wb, file.path(ds_out, "DEG_TP53mt_vs_TP53wt.xlsx"), overwrite = TRUE)

        ds_summary$TP53mt_vs_TP53wt_total <- nrow(deg1)
        ds_summary$TP53mt_vs_TP53wt_sig <- sum(deg1$adj.P.Val < 0.05, na.rm = TRUE)
        ds_summary$TP53mt_vs_TP53wt_up <- sum(deg1$adj.P.Val < 0.05 & deg1$logFC > 0, na.rm = TRUE)
        ds_summary$TP53mt_vs_TP53wt_down <- sum(deg1$adj.P.Val < 0.05 & deg1$logFC < 0, na.rm = TRUE)
    } else {
        ds_summary$TP53mt_vs_TP53wt_total <- NA
        ds_summary$TP53mt_vs_TP53wt_sig <- NA
        ds_summary$TP53mt_vs_TP53wt_up <- NA
        ds_summary$TP53mt_vs_TP53wt_down <- NA
    }

    # ========================================================================
    # Comparison 2: MUT_GOF vs MUT_LOF
    # ========================================================================
    cat("\n  --- Comparison 2: MUT_GOF vs MUT_LOF ---\n")

    tp53_mut <- tp53_ds %>% filter(mt == 1)
    gof_samples <- tp53_mut$sample_id[tp53_mut$GOF == 1]
    lof_samples <- tp53_mut$sample_id[tp53_mut$LOF == 1]

    cat("    GOF samples:", length(gof_samples), " LOF samples:", length(lof_samples), "\n")

    group2_ids <- c(gof_samples, lof_samples)
    group2 <- setNames(
        c(rep("MUT_GOF", length(gof_samples)), rep("MUT_LOF", length(lof_samples))),
        group2_ids
    )
    group2 <- factor(group2, levels = c("MUT_LOF", "MUT_GOF"))

    deg2 <- run_limma_group_comparison(mat, group2, batch_fac, purity_vec, sa_df)

    if (!is.null(deg2)) {
        fwrite(deg2, file.path(ds_out, "DEG_MUT_GOF_vs_MUT_LOF.csv"))

        wb <- createWorkbook()
        addWorksheet(wb, "GOF_vs_LOF")
        writeData(wb, 1, deg2)
        saveWorkbook(wb, file.path(ds_out, "DEG_MUT_GOF_vs_MUT_LOF.xlsx"), overwrite = TRUE)

        ds_summary$GOF_vs_LOF_total <- nrow(deg2)
        ds_summary$GOF_vs_LOF_sig <- sum(deg2$adj.P.Val < 0.05, na.rm = TRUE)
        ds_summary$GOF_vs_LOF_up <- sum(deg2$adj.P.Val < 0.05 & deg2$logFC > 0, na.rm = TRUE)
        ds_summary$GOF_vs_LOF_down <- sum(deg2$adj.P.Val < 0.05 & deg2$logFC < 0, na.rm = TRUE)
    } else {
        ds_summary$GOF_vs_LOF_total <- NA
        ds_summary$GOF_vs_LOF_sig <- NA
        ds_summary$GOF_vs_LOF_up <- NA
        ds_summary$GOF_vs_LOF_down <- NA
    }

    # ========================================================================
    # Comparison 3: Hotspots vs MUT_LOF
    # ========================================================================
    cat("\n  --- Comparison 3: Hotspots vs MUT_LOF ---\n")

    hotspot_samples <- tp53_mut$sample_id[tp53_mut$hotspot == 1]

    cat("    Hotspot samples:", length(hotspot_samples), " LOF samples:", length(lof_samples), "\n")

    group3_ids <- c(hotspot_samples, lof_samples)
    # Remove duplicates (a hotspot might also be LOF classified)
    group3_labels <- setNames(
        c(rep("Hotspot", length(hotspot_samples)), rep("MUT_LOF", length(lof_samples))),
        group3_ids
    )
    # If a sample appears in both, keep Hotspot classification
    group3_labels <- group3_labels[!duplicated(names(group3_labels))]
    group3 <- factor(group3_labels, levels = c("MUT_LOF", "Hotspot"))

    deg3 <- run_limma_group_comparison(mat, group3, batch_fac, purity_vec, sa_df)

    if (!is.null(deg3)) {
        fwrite(deg3, file.path(ds_out, "DEG_Hotspot_vs_MUT_LOF.csv"))

        wb <- createWorkbook()
        addWorksheet(wb, "Hotspot_vs_LOF")
        writeData(wb, 1, deg3)
        saveWorkbook(wb, file.path(ds_out, "DEG_Hotspot_vs_MUT_LOF.xlsx"), overwrite = TRUE)

        ds_summary$Hotspot_vs_LOF_total <- nrow(deg3)
        ds_summary$Hotspot_vs_LOF_sig <- sum(deg3$adj.P.Val < 0.05, na.rm = TRUE)
        ds_summary$Hotspot_vs_LOF_up <- sum(deg3$adj.P.Val < 0.05 & deg3$logFC > 0, na.rm = TRUE)
        ds_summary$Hotspot_vs_LOF_down <- sum(deg3$adj.P.Val < 0.05 & deg3$logFC < 0, na.rm = TRUE)
    } else {
        ds_summary$Hotspot_vs_LOF_total <- NA
        ds_summary$Hotspot_vs_LOF_sig <- NA
        ds_summary$Hotspot_vs_LOF_up <- NA
        ds_summary$Hotspot_vs_LOF_down <- NA
    }

    # ========================================================================
    # Comparison 4: MUT_GOF vs TP53wt
    # ========================================================================
    cat("\n  --- Comparison 4: MUT_GOF vs TP53wt ---\n")

    wt_samples <- tp53_ds$sample_id[tp53_ds$wt == 1]
    cat("    GOF samples:", length(gof_samples), " TP53wt samples:", length(wt_samples), "\n")

    group4_ids <- c(gof_samples, wt_samples)
    group4 <- setNames(
        c(rep("MUT_GOF", length(gof_samples)), rep("TP53wt", length(wt_samples))),
        group4_ids
    )
    group4 <- factor(group4, levels = c("TP53wt", "MUT_GOF"))

    deg4 <- run_limma_group_comparison(mat, group4, batch_fac, purity_vec, sa_df)

    if (!is.null(deg4)) {
        fwrite(deg4, file.path(ds_out, "DEG_MUT_GOF_vs_TP53wt.csv"))
        wb <- createWorkbook()
        addWorksheet(wb, "GOF_vs_TP53wt")
        writeData(wb, 1, deg4)
        saveWorkbook(wb, file.path(ds_out, "DEG_MUT_GOF_vs_TP53wt.xlsx"), overwrite = TRUE)
        ds_summary$GOF_vs_TP53wt_total <- nrow(deg4)
        ds_summary$GOF_vs_TP53wt_sig <- sum(deg4$adj.P.Val < 0.05, na.rm = TRUE)
        ds_summary$GOF_vs_TP53wt_up <- sum(deg4$adj.P.Val < 0.05 & deg4$logFC > 0, na.rm = TRUE)
        ds_summary$GOF_vs_TP53wt_down <- sum(deg4$adj.P.Val < 0.05 & deg4$logFC < 0, na.rm = TRUE)
    } else {
        ds_summary$GOF_vs_TP53wt_total <- NA
        ds_summary$GOF_vs_TP53wt_sig <- NA
        ds_summary$GOF_vs_TP53wt_up <- NA
        ds_summary$GOF_vs_TP53wt_down <- NA
    }

    # ========================================================================
    # Comparison 5: MUT_LOF vs TP53wt
    # ========================================================================
    cat("\n  --- Comparison 5: MUT_LOF vs TP53wt ---\n")

    cat("    LOF samples:", length(lof_samples), " TP53wt samples:", length(wt_samples), "\n")

    group5_ids <- c(lof_samples, wt_samples)
    group5 <- setNames(
        c(rep("MUT_LOF", length(lof_samples)), rep("TP53wt", length(wt_samples))),
        group5_ids
    )
    group5 <- factor(group5, levels = c("TP53wt", "MUT_LOF"))

    deg5 <- run_limma_group_comparison(mat, group5, batch_fac, purity_vec, sa_df)

    if (!is.null(deg5)) {
        fwrite(deg5, file.path(ds_out, "DEG_MUT_LOF_vs_TP53wt.csv"))
        wb <- createWorkbook()
        addWorksheet(wb, "LOF_vs_TP53wt")
        writeData(wb, 1, deg5)
        saveWorkbook(wb, file.path(ds_out, "DEG_MUT_LOF_vs_TP53wt.xlsx"), overwrite = TRUE)
        ds_summary$LOF_vs_TP53wt_total <- nrow(deg5)
        ds_summary$LOF_vs_TP53wt_sig <- sum(deg5$adj.P.Val < 0.05, na.rm = TRUE)
        ds_summary$LOF_vs_TP53wt_up <- sum(deg5$adj.P.Val < 0.05 & deg5$logFC > 0, na.rm = TRUE)
        ds_summary$LOF_vs_TP53wt_down <- sum(deg5$adj.P.Val < 0.05 & deg5$logFC < 0, na.rm = TRUE)
    } else {
        ds_summary$LOF_vs_TP53wt_total <- NA
        ds_summary$LOF_vs_TP53wt_sig <- NA
        ds_summary$LOF_vs_TP53wt_up <- NA
        ds_summary$LOF_vs_TP53wt_down <- NA
    }

    # ========================================================================
    # Comparison 6: Hotspots vs TP53wt
    # ========================================================================
    cat("\n  --- Comparison 6: Hotspots vs TP53wt ---\n")

    cat("    Hotspot samples:", length(hotspot_samples), " TP53wt samples:", length(wt_samples), "\n")

    group6_ids <- c(hotspot_samples, wt_samples)
    group6 <- setNames(
        c(rep("Hotspot", length(hotspot_samples)), rep("TP53wt", length(wt_samples))),
        group6_ids
    )
    group6 <- factor(group6, levels = c("TP53wt", "Hotspot"))

    deg6 <- run_limma_group_comparison(mat, group6, batch_fac, purity_vec, sa_df)

    if (!is.null(deg6)) {
        fwrite(deg6, file.path(ds_out, "DEG_Hotspot_vs_TP53wt.csv"))
        wb <- createWorkbook()
        addWorksheet(wb, "Hotspot_vs_TP53wt")
        writeData(wb, 1, deg6)
        saveWorkbook(wb, file.path(ds_out, "DEG_Hotspot_vs_TP53wt.xlsx"), overwrite = TRUE)
        ds_summary$Hotspot_vs_TP53wt_total <- nrow(deg6)
        ds_summary$Hotspot_vs_TP53wt_sig <- sum(deg6$adj.P.Val < 0.05, na.rm = TRUE)
        ds_summary$Hotspot_vs_TP53wt_up <- sum(deg6$adj.P.Val < 0.05 & deg6$logFC > 0, na.rm = TRUE)
        ds_summary$Hotspot_vs_TP53wt_down <- sum(deg6$adj.P.Val < 0.05 & deg6$logFC < 0, na.rm = TRUE)
    } else {
        ds_summary$Hotspot_vs_TP53wt_total <- NA
        ds_summary$Hotspot_vs_TP53wt_sig <- NA
        ds_summary$Hotspot_vs_TP53wt_up <- NA
        ds_summary$Hotspot_vs_TP53wt_down <- NA
    }

    # ========================================================================
    # Comparison 7: DN vs TP53wt
    # ========================================================================
    cat("\n  --- Comparison 7: DN vs TP53wt ---\n")

    dn_samples <- tp53_mut$sample_id[tp53_mut$DN == 1]
    cat("    DN samples:", length(dn_samples), " TP53wt samples:", length(wt_samples), "\n")

    group7_ids <- c(dn_samples, wt_samples)
    group7 <- setNames(
        c(rep("DN", length(dn_samples)), rep("TP53wt", length(wt_samples))),
        group7_ids
    )
    group7 <- factor(group7, levels = c("TP53wt", "DN"))

    deg7 <- run_limma_group_comparison(mat, group7, batch_fac, purity_vec, sa_df)

    if (!is.null(deg7)) {
        fwrite(deg7, file.path(ds_out, "DEG_DN_vs_TP53wt.csv"))
        wb <- createWorkbook()
        addWorksheet(wb, "DN_vs_TP53wt")
        writeData(wb, 1, deg7)
        saveWorkbook(wb, file.path(ds_out, "DEG_DN_vs_TP53wt.xlsx"), overwrite = TRUE)
        ds_summary$DN_vs_TP53wt_total <- nrow(deg7)
        ds_summary$DN_vs_TP53wt_sig <- sum(deg7$adj.P.Val < 0.05, na.rm = TRUE)
        ds_summary$DN_vs_TP53wt_up <- sum(deg7$adj.P.Val < 0.05 & deg7$logFC > 0, na.rm = TRUE)
        ds_summary$DN_vs_TP53wt_down <- sum(deg7$adj.P.Val < 0.05 & deg7$logFC < 0, na.rm = TRUE)
    } else {
        ds_summary$DN_vs_TP53wt_total <- NA
        ds_summary$DN_vs_TP53wt_sig <- NA
        ds_summary$DN_vs_TP53wt_up <- NA
        ds_summary$DN_vs_TP53wt_down <- NA
    }

    # ========================================================================
    # Comparison 8: Non-DN vs TP53wt
    # ========================================================================
    cat("\n  --- Comparison 8: Non-DN vs TP53wt ---\n")

    nondn_samples <- tp53_mut$sample_id[tp53_mut$non_DN == 1]
    cat("    Non-DN samples:", length(nondn_samples), " TP53wt samples:", length(wt_samples), "\n")

    group8_ids <- c(nondn_samples, wt_samples)
    group8 <- setNames(
        c(rep("NonDN", length(nondn_samples)), rep("TP53wt", length(wt_samples))),
        group8_ids
    )
    group8 <- factor(group8, levels = c("TP53wt", "NonDN"))

    deg8 <- run_limma_group_comparison(mat, group8, batch_fac, purity_vec, sa_df)

    if (!is.null(deg8)) {
        fwrite(deg8, file.path(ds_out, "DEG_NonDN_vs_TP53wt.csv"))
        wb <- createWorkbook()
        addWorksheet(wb, "NonDN_vs_TP53wt")
        writeData(wb, 1, deg8)
        saveWorkbook(wb, file.path(ds_out, "DEG_NonDN_vs_TP53wt.xlsx"), overwrite = TRUE)
        ds_summary$NonDN_vs_TP53wt_total <- nrow(deg8)
        ds_summary$NonDN_vs_TP53wt_sig <- sum(deg8$adj.P.Val < 0.05, na.rm = TRUE)
        ds_summary$NonDN_vs_TP53wt_up <- sum(deg8$adj.P.Val < 0.05 & deg8$logFC > 0, na.rm = TRUE)
        ds_summary$NonDN_vs_TP53wt_down <- sum(deg8$adj.P.Val < 0.05 & deg8$logFC < 0, na.rm = TRUE)
    } else {
        ds_summary$NonDN_vs_TP53wt_total <- NA
        ds_summary$NonDN_vs_TP53wt_sig <- NA
        ds_summary$NonDN_vs_TP53wt_up <- NA
        ds_summary$NonDN_vs_TP53wt_down <- NA
    }

    # ========================================================================
    # Comparison 9: DN vs non-DN
    # ========================================================================
    cat("\n  --- Comparison 9: DN vs non-DN ---\n")

    cat("    DN samples:", length(dn_samples), " Non-DN samples:", length(nondn_samples), "\n")

    group9_ids <- c(dn_samples, nondn_samples)
    group9_labels <- setNames(
        c(rep("DN", length(dn_samples)), rep("NonDN", length(nondn_samples))),
        group9_ids
    )
    # Handle possible overlap: a sample could be both DN and non-DN (shouldn't happen, but safe)
    group9_labels <- group9_labels[!duplicated(names(group9_labels))]
    group9 <- factor(group9_labels, levels = c("NonDN", "DN"))

    deg9 <- run_limma_group_comparison(mat, group9, batch_fac, purity_vec, sa_df)

    if (!is.null(deg9)) {
        fwrite(deg9, file.path(ds_out, "DEG_DN_vs_NonDN.csv"))
        wb <- createWorkbook()
        addWorksheet(wb, "DN_vs_NonDN")
        writeData(wb, 1, deg9)
        saveWorkbook(wb, file.path(ds_out, "DEG_DN_vs_NonDN.xlsx"), overwrite = TRUE)
        ds_summary$DN_vs_NonDN_total <- nrow(deg9)
        ds_summary$DN_vs_NonDN_sig <- sum(deg9$adj.P.Val < 0.05, na.rm = TRUE)
        ds_summary$DN_vs_NonDN_up <- sum(deg9$adj.P.Val < 0.05 & deg9$logFC > 0, na.rm = TRUE)
        ds_summary$DN_vs_NonDN_down <- sum(deg9$adj.P.Val < 0.05 & deg9$logFC < 0, na.rm = TRUE)
    } else {
        ds_summary$DN_vs_NonDN_total <- NA
        ds_summary$DN_vs_NonDN_sig <- NA
        ds_summary$DN_vs_NonDN_up <- NA
        ds_summary$DN_vs_NonDN_down <- NA
    }

    all_deg_summary[[ds_folder]] <- ds_summary
    cat("\n")
}

# ==============================================================================
# Section 5: Summary Statistics
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
cat("\n=== Differential Gene Summary (padj < 0.05) ===\n\n")
print(as.data.frame(summary_df), row.names = FALSE)

cat("\n====================================================================\n")
cat("Pipeline Complete!\n")
cat("All results saved to:", output_dir, "\n")
cat("====================================================================\n")
