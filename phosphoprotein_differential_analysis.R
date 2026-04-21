################################################################################
# TP53 Phosphoprotein Differential Phosphosite Analysis (limma)
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
# Protein normalization applied (phospho - protein)
#
# Methodology adapted from phosphoprotein_DPS_analysis/ pipeline
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
    library(vroom)
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

output_dir <- file.path(base_path, "phosphoprotein_differential_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

min_frac_complete <- 0.50  # phosphosite data is sparser
min_per_group <- 5

cat("====================================================================\n")
cat("TP53 Phosphoprotein Differential Phosphosite Analysis\n")
cat("====================================================================\n\n")

# ==============================================================================
# Section 2: Utility Functions (adapted from phosphoprotein_DPS_analysis)
# ==============================================================================

#' Standardize site ID to GENE_S123 / GENE_T123 / GENE_Y123 format
.std_site_id <- function(x) {
    x <- toupper(as.character(x))
    x <- gsub("[:\\s]", "_", x)
    x <- gsub("_([STY])(\\d+)[A-Z]*$", "_\\1\\2", x)
    x
}

#' Infer site ID from phosphoprotein metadata columns
.infer_site_id <- function(df) {
    g <- if ("GENE_SYMBOL" %in% names(df)) as.character(df$GENE_SYMBOL) else rep(NA_character_, nrow(df))
    nm <- if ("NAME" %in% names(df)) as.character(df$NAME) else rep(NA_character_, nrow(df))
    p1 <- if ("PHOSPHOSITES" %in% names(df)) as.character(df$PHOSPHOSITES) else rep(NA_character_, nrow(df))
    p2 <- if ("PHOSPHOSITE" %in% names(df)) as.character(df$PHOSPHOSITE) else rep(NA_character_, nrow(df))
    ent <- if ("ENTITY_STABLE_ID" %in% names(df)) as.character(df$ENTITY_STABLE_ID) else rep(NA_character_, nrow(df))

    # 1) from NAME
    id_from_name <- ifelse(!is.na(nm) & grepl("_[STY]\\d+", nm, ignore.case = TRUE),
        sub(":.*$", "", nm), NA_character_
    )

    site_token <- function(x) {
        y <- NA_character_
        if (!is.na(x) && nzchar(x)) {
            m <- regmatches(x, regexpr("([STY])(\\d+)", x, ignore.case = TRUE))
            if (length(m) == 1 && nzchar(m)) y <- toupper(m)
        }
        y
    }

    # 2) from GENE + PHOSPHOSITES
    ps <- vapply(seq_along(p1), function(i) {
        s <- if (!is.na(p1[i]) && nzchar(p1[i])) p1[i] else p2[i]
        site_token(s)
    }, character(1))
    id_from_gp <- ifelse(!is.na(g) & nzchar(g) & !is.na(ps) & nzchar(ps),
        paste0(g, "_", ps), NA_character_
    )

    # 3) from ENTITY_STABLE_ID
    ent_base <- toupper(sub(":.*$", "", ent))
    ent_gene <- ifelse(!is.na(g) & nzchar(g), toupper(g), sub("_.*$", "", ent_base))
    ent_site <- vapply(ent_base, function(x) {
        m <- regmatches(x, regexpr("([STY])(\\d+)", x, ignore.case = TRUE))
        if (length(m) == 1 && nzchar(m)) toupper(m) else NA_character_
    }, character(1))
    id_from_ent <- ifelse(!is.na(ent_gene) & nzchar(ent_gene) & !is.na(ent_site) & nzchar(ent_site),
        paste0(ent_gene, "_", ent_site), NA_character_
    )

    cand <- id_from_name
    need <- is.na(cand) | !nzchar(cand)
    cand[need] <- id_from_gp[need]
    need <- is.na(cand) | !nzchar(cand)
    cand[need] <- id_from_ent[need]

    .std_site_id(cand)
}

#' Load site-level phosphosite matrix with protein normalization
load_phosphosite_matrix <- function(ds_dir, protein_adjust = TRUE) {
    fp <- file.path(ds_dir, "data_phosphoprotein_quantification.txt")
    if (!file.exists(fp)) stop("Cannot find phosphoprotein file: ", fp)
    cat("  Reading phosphosite matrix:", basename(fp), "\n")

    df <- suppressMessages(vroom::vroom(fp, delim = "\t", col_types = vroom::cols(
        .default = "d",
        ENTITY_STABLE_ID = "c", NAME = "c", DESCRIPTION = "c",
        GENE_SYMBOL = "c", PHOSPHOSITES = "c", PHOSPHOSITE = "c"
    )))
    df <- as.data.frame(df, check.names = FALSE)

    # Detect annotation vs sample columns
    anno_cols <- intersect(
        c("ENTITY_STABLE_ID", "NAME", "DESCRIPTION", "GENE_SYMBOL", "PHOSPHOSITES", "PHOSPHOSITE"),
        colnames(df)
    )
    samp_cols <- setdiff(colnames(df), anno_cols)
    stopifnot(length(samp_cols) > 0)

    # Generate site IDs
    site_id <- .infer_site_id(df)
    keep <- !is.na(site_id) & nzchar(site_id)
    df_keep <- df[keep, , drop = FALSE]
    site_id <- site_id[keep]

    M_site <- as.matrix(df_keep[, samp_cols, drop = FALSE])
    storage.mode(M_site) <- "numeric"
    
    # Aggregate duplicate site IDs (e.g. multiple peptides for the same site)
    if (anyDuplicated(site_id)) {
        split_idx <- split(seq_len(nrow(M_site)), site_id)
        M_site_agg <- do.call(rbind, lapply(split_idx, function(idx) {
            if (length(idx) == 1) return(M_site[idx, ])
            apply(M_site[idx, , drop = FALSE], 2, function(x) stats::median(x, na.rm = TRUE))
        }))
        rownames(M_site_agg) <- names(split_idx)
        site_id <- rownames(M_site_agg)
        M_site <- M_site_agg
    } else {
        rownames(M_site) <- site_id
    }
    colnames(M_site) <- samp_cols

    # Protein normalization: subtract protein abundance from phosphosite
    if (isTRUE(protein_adjust)) {
        prot_fp <- file.path(ds_dir, "data_protein_quantification.txt")
        if (!file.exists(prot_fp)) {
            cat("  WARNING: protein file not found, skipping normalization\n")
        } else {
            cat("  Applying protein normalization\n")
            prot_df <- suppressMessages(vroom::vroom(prot_fp, delim = "\t"))
            prot_df <- as.data.frame(prot_df, check.names = FALSE)

            gcol <- which(tolower(colnames(prot_df)) %in%
                tolower(c("Composite.Element.REF", "Gene", "GENE", "GENE_SYMBOL",
                           "GENE_NAME", "Hugo_Symbol")))[1]
            if (!is.na(gcol)) {
                prot_genes <- as.character(prot_df[[gcol]])
                prot_mat <- as.matrix(prot_df[, setdiff(colnames(prot_df), colnames(prot_df)[gcol]), drop = FALSE])
                rownames(prot_mat) <- prot_genes
                storage.mode(prot_mat) <- "numeric"

                gene_of_site <- as.character(df_keep$GENE_SYMBOL)

                common <- intersect(colnames(M_site), colnames(prot_mat))
                if (length(common) >= 2) {
                    M_site <- M_site[, common, drop = FALSE]
                    prot_mat <- prot_mat[, common, drop = FALSE]

                    idx <- match(gene_of_site, rownames(prot_mat))
                    hasp <- !is.na(idx)
                    if (any(hasp)) {
                        M_site[hasp, ] <- M_site[hasp, , drop = FALSE] - prot_mat[idx[hasp], , drop = FALSE]
                        cat("    Normalized", sum(hasp), "of", nrow(M_site), "phosphosites\n")
                    }
                }
            }
        }
    }

    cat("  Phosphosite matrix:", nrow(M_site), "sites x", ncol(M_site), "samples\n")
    M_site
}

#' Impute and filter phosphosite matrix (gene-median imputation)
impute_and_filter_phospho <- function(mat, min_frac = 0.50) {
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

#' Get batch factor for phosphoprotein data
get_batch_factor_phospho <- function(ds_dir, sample_ids) {
    meta_fp <- file.path(ds_dir, "data_clinical_sample.txt")
    if (file.exists(meta_fp)) {
        meta <- suppressMessages(
            readr::read_tsv(meta_fp, show_col_types = FALSE, comment = "#")
        ) |> as.data.frame()
        id_cols <- intersect(c("SAMPLE_ID", "sample_id", "Sample_ID"), names(meta))
        if (length(id_cols)) {
            meta <- dplyr::rename(meta, SAMPLE_ID = !!id_cols[1])
            meta$SAMPLE_ID <- as.character(meta$SAMPLE_ID)
            meta <- meta[match(sample_ids, meta$SAMPLE_ID), , drop = FALSE]
            rownames(meta) <- sample_ids
            cand <- c("TMT_PLEX", "EXPERIMENT")
            hit <- intersect(cand, colnames(meta))
            for (cn in hit) {
                fac <- sanitize_batch_levels(meta[[cn]])
                if (nlevels(fac) >= 2 && sum(!is.na(fac)) >= 3) {
                    names(fac) <- sample_ids
                    return(list(name = cn, fac = fac))
                }
            }
        }
    }
    # Fallback: TMT_phos.csv
    tmt_fp <- file.path(ds_dir, "TMT_phos.csv")
    if (file.exists(tmt_fp)) {
        tmt <- suppressMessages(readr::read_csv(tmt_fp, show_col_types = FALSE)) |> as.data.frame()
        cn <- trimws(names(tmt))
        names(tmt) <- cn
        run_hits <- grep("^Run\\s*Metadata\\s*ID$", cn, ignore.case = TRUE, value = TRUE)
        tmt_cols <- grep("^tmt_", cn, ignore.case = TRUE, value = TRUE)
        if (length(run_hits) && length(tmt_cols)) {
            run_col <- run_hits[1]
            plex_by_sample <- list()
            for (i in seq_len(nrow(tmt))) {
                run_id <- as.character(tmt[[run_col]][i])
                if (!nzchar(run_id) || is.na(run_id)) next
                plex2 <- sub("^\\s*(\\d{2}).*$", "\\1", run_id)
                for (tc in tmt_cols) {
                    cell <- as.character(tmt[[tc]][i])
                    if (is.na(cell) || !nzchar(cell)) next
                    sid <- trimws(sub("\\r?\\n.*$", "", cell))
                    if (nzchar(sid) && is.null(plex_by_sample[[sid]])) plex_by_sample[[sid]] <- plex2
                }
            }
            if (length(plex_by_sample)) {
                v <- setNames(rep(NA_character_, length(sample_ids)), sample_ids)
                mm <- intersect(names(plex_by_sample), sample_ids)
                if (length(mm)) v[mm] <- unlist(plex_by_sample[mm], use.names = FALSE)
                fac2 <- sanitize_batch_levels(v)
                if (nlevels(fac2) >= 2 && sum(!is.na(fac2)) >= 3) {
                    names(fac2) <- sample_ids
                    return(list(name = "TMT_phos.csv", fac = fac2))
                }
            }
        }
    }
    NULL
}

#' Helper functions for purity extraction
.norm_names <- function(x) toupper(gsub("[^A-Za-z0-9]+", "_", x))

.to01 <- function(v) {
    v <- suppressWarnings(as.numeric(v))
    if (sum(is.finite(v) & v > 1, na.rm = TRUE) > sum(is.finite(v) & v <= 1, na.rm = TRUE)) v <- v / 100
    pmin(pmax(v, 0), 1)
}

.median_from_semicolon <- function(x_chr) {
    vv <- suppressWarnings(as.numeric(unlist(strsplit(as.character(x_chr), ";"))))
    vv <- vv[is.finite(vv)]
    if (!length(vv)) return(NA_real_)
    stats::median(vv)
}

.z_no_impute <- function(v) {
    v <- suppressWarnings(as.numeric(v))
    fin <- is.finite(v)
    mu <- mean(v[fin], na.rm = TRUE)
    sdv <- stats::sd(v[fin], na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    out <- (v - mu) / sdv
    out[!fin] <- NA_real_
    out
}

#' Get tumor purity (dataset-specific logic from phosphoprotein_DPS_analysis)
get_purity_covariate <- function(ds_id, ds_dir, sample_ids) {
    samp_fp <- file.path(ds_dir, "data_clinical_sample.txt")
    pat_fp <- file.path(ds_dir, "data_clinical_patient.txt")
    samp <- if (file.exists(samp_fp)) suppressMessages(readr::read_tsv(samp_fp, show_col_types = FALSE, comment = "#")) else NULL
    if (!is.null(samp)) { samp <- as.data.frame(samp); names(samp) <- .norm_names(names(samp)) }
    purity <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)

    if (ds_id == "brca_cptac_2020") {
        if (file.exists(pat_fp) && !is.null(samp)) {
            pat <- suppressMessages(readr::read_tsv(pat_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
            names(pat) <- .norm_names(names(pat))
            sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
            pid_in_samp <- intersect(c("PATIENT_ID", "PATIENT"), names(samp))[1]
            pid_in_pat <- intersect(c("PATIENT_ID", "PATIENT"), names(pat))[1]
            if (!is.na(sid) && !is.na(pid_in_samp) && !is.na(pid_in_pat) && "ESTIMATE_TUMORPURITY" %in% names(pat)) {
                map_pt <- setNames(as.character(samp[[pid_in_samp]]), samp[[sid]])
                pt_pur <- setNames(.to01(pat[["ESTIMATE_TUMORPURITY"]]), as.character(pat[[pid_in_pat]]))
                purity[] <- unname(pt_pur[map_pt[sample_ids]])
            }
        }
    } else if (ds_id == "luad_cptac_2020") {
        if (!is.null(samp)) {
            sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
            if (!is.na(sid) && "TUMOR_PURITY_BYESTIMATE_RNASEQ" %in% names(samp))
                purity[] <- .to01(samp[["TUMOR_PURITY_BYESTIMATE_RNASEQ"]][match(sample_ids, samp[[sid]])])
        }
    } else if (ds_id == "paad_cptac_2021") {
        if (!is.null(samp)) {
            sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
            if (!is.na(sid) && "NEOPLASTIC_CELLULARITY" %in% names(samp)) {
                medv <- vapply(as.character(samp[["NEOPLASTIC_CELLULARITY"]]), .median_from_semicolon, numeric(1))
                purity[] <- .to01(medv[match(sample_ids, samp[[sid]])])
            }
        }
    } else if (ds_id == "ucec_cptac_2020") {
        if (!is.null(samp)) {
            sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
            if (!is.na(sid) && "PURITY_CANCER" %in% names(samp)) {
                purity[] <- .to01(suppressWarnings(as.numeric(samp[["PURITY_CANCER"]]))[match(sample_ids, samp[[sid]])])
            }
        }
    } else {
        # Generic fallback: try TUMOR_PURITY from patient file
        if (file.exists(pat_fp)) {
            pat <- tryCatch(readr::read_tsv(pat_fp, comment = "#", show_col_types = FALSE), error = function(e) NULL)
            if (!is.null(pat)) {
                pat <- as.data.frame(pat); names(pat) <- .norm_names(names(pat))
                purity_cols <- c("TUMOR_PURITY", "PURITY")
                hit <- intersect(purity_cols, names(pat))
                if (length(hit) && !is.null(samp)) {
                    sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
                    pid_in_samp <- intersect(c("PATIENT_ID", "PATIENT"), names(samp))[1]
                    pid_in_pat <- intersect(c("PATIENT_ID", "PATIENT"), names(pat))[1]
                    if (!is.na(sid) && !is.na(pid_in_samp) && !is.na(pid_in_pat)) {
                        map_pt <- setNames(as.character(samp[[pid_in_samp]]), samp[[sid]])
                        pt_pur <- setNames(suppressWarnings(as.numeric(pat[[hit[1]]])), as.character(pat[[pid_in_pat]]))
                        purity[] <- unname(pt_pur[map_pt[sample_ids]])
                    }
                }
            }
        }
    }
    purity
}

#' Get sex and age covariates (from patient file, following phosphoprotein pipeline)
get_sex_age_covariates <- function(ds_dir, sample_ids) {
    samp_fp <- file.path(ds_dir, "data_clinical_sample.txt")
    pat_fp <- file.path(ds_dir, "data_clinical_patient.txt")
    empty_df <- data.frame(sex = rep(NA_real_, length(sample_ids)),
                           age = rep(NA_real_, length(sample_ids)),
                           row.names = sample_ids, stringsAsFactors = FALSE)
    if (!file.exists(samp_fp) || !file.exists(pat_fp)) return(empty_df)

    samp <- suppressMessages(readr::read_tsv(samp_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
    pat <- suppressMessages(readr::read_tsv(pat_fp, show_col_types = FALSE, comment = "#")) |> as.data.frame()
    names(samp) <- .norm_names(names(samp))
    names(pat) <- .norm_names(names(pat))

    sid <- intersect(c("SAMPLE_ID", "SAMPLE"), names(samp))[1]
    pid_samp <- intersect(c("PATIENT_ID", "PATIENT"), names(samp))[1]
    pid_pat <- intersect(c("PATIENT_ID", "PATIENT"), names(pat))[1]
    if (is.na(sid) || is.na(pid_samp) || is.na(pid_pat)) return(empty_df)
    map_pt <- setNames(as.character(samp[[pid_samp]]), samp[[sid]])

    # SEX: Male=1, Female=0
    sex_col <- intersect(c("SEX", "GENDER"), names(pat))[1]
    sex <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
    if (!is.na(sex_col)) {
        raw <- toupper(as.character(pat[[sex_col]]))
        val <- ifelse(grepl("^M", raw), 1, ifelse(grepl("^F", raw), 0, NA_real_))
        names(val) <- as.character(pat[[pid_pat]])
        sex[] <- unname(val[map_pt[sample_ids]])
    }

    # AGE (z-scored, no imputation)
    age_col <- intersect(c("AGE", "AGE_AT_DIAGNOSIS", "AGE_AT_INDEX"), names(pat))[1]
    age <- setNames(rep(NA_real_, length(sample_ids)), sample_ids)
    if (!is.na(age_col)) {
        v <- suppressWarnings(as.numeric(pat[[age_col]]))
        names(v) <- as.character(pat[[pid_pat]])
        age_raw <- unname(v[map_pt[sample_ids]])
        age[] <- .z_no_impute(age_raw)
    }

    data.frame(sex = as.numeric(sex), age = as.numeric(age),
               row.names = sample_ids, stringsAsFactors = FALSE)
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
# Section 3: Limma DEG Function for Group Comparison (Phosphosite)
# ==============================================================================

run_limma_group_comparison <- function(mat, group_vec,
                                       batch_fac = NULL,
                                       purity_vec = NULL,
                                       sa_df = NULL) {
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
            if (sum(is.finite(s)) >= length(so) * 0.8) {
                DF$sex <- s
            }
        }
        if ("age" %in% colnames(sa_df)) {
            a <- suppressWarnings(as.numeric(sa_df[so, "age"]))
            if (sum(is.finite(a)) >= length(so) * 0.8) {
                DF$age <- a
            }
        }
    }

    # Clean
    all_na <- vapply(DF, function(z) all(is.na(z)), logical(1))
    if (any(all_na)) DF <- DF[, !all_na, drop = FALSE]
    DF <- coerce_covariates_safely(DF)

    if (!("group" %in% colnames(DF))) {
        cat("    [ERROR] Group variable dropped\n")
        return(NULL)
    }

    # Remove rows with any NA (critical: model.matrix silently drops NA rows)
    ok <- complete.cases(DF)
    DF <- DF[ok, , drop = FALSE]

    # Re-check group sizes after NA removal
    if (nrow(DF) < 2 * min_per_group) {
        cat("    [SKIP] Too few complete-case samples:", nrow(DF), "\n")
        return(NULL)
    }
    tab2 <- table(DF$group)
    if (any(tab2 < min_per_group)) {
        cat("    [SKIP] Group too small after NA removal:", paste(names(tab2), "=", tab2, collapse = ", "), "\n")
        return(NULL)
    }

    M <- M[, rownames(DF), drop = FALSE]

    des <- model.matrix(~ 0 + ., data = DF)
    lvls <- levels(DF$group)
    group_cols <- paste0("group", lvls)
    if (!all(group_cols %in% colnames(des))) {
        cat("    [ERROR] Group columns not found in design matrix\n")
        return(NULL)
    }

    # Saturation protection: drop covariates if residual df <= 0
    rnk <- qr(des)$rank
    res_df <- nrow(des) - rnk
    if (res_df <= 0) {
        drop_order <- intersect(c("batch", "sex", "age", "purity"), colnames(DF))
        for (cv in drop_order) {
            DF[[cv]] <- NULL
            cat("    [saturation] dropping covariate:", cv, "\n")
            des <- model.matrix(~ 0 + ., data = DF)
            rnk <- qr(des)$rank
            res_df <- nrow(des) - rnk
            if (res_df > 0) break
        }
        if (res_df <= 0) {
            cat("    [SKIP] No residual df after removing covariates\n")
            return(NULL)
        }
        # Re-derive group columns after design change
        lvls <- levels(DF$group)
        group_cols <- paste0("group", lvls)
        M <- M[, rownames(DF), drop = FALSE]
    }

    contrast_str <- paste0(group_cols[2], " - ", group_cols[1])
    contrast_mat <- limma::makeContrasts(contrasts = contrast_str, levels = des)

    nobs <- rowSums(is.finite(M))
    keep_rows <- nobs >= (rnk + 1L)
    if (sum(keep_rows) == 0) {
        cat("    [SKIP] No phosphosites pass df filter\n")
        return(NULL)
    }
    M <- M[keep_rows, , drop = FALSE]

    fit <- limma::lmFit(M, design = des)
    fit2 <- limma::contrasts.fit(fit, contrast_mat)
    eb <- limma::eBayes(fit2, trend = TRUE)

    tbl <- limma::topTable(eb, number = Inf, sort.by = "P")
    tbl$phosphosite <- rownames(tbl)
    tbl <- tbl[, c("phosphosite", "logFC", "t", "P.Value", "adj.P.Val", "B")]

    cat("    Results:", nrow(tbl), "phosphosites tested,",
        sum(tbl$adj.P.Val < 0.05, na.rm = TRUE), "significant (padj < 0.05)\n")
    tbl
}

# ==============================================================================
# Section 4: Helper to run one comparison and save results
# ==============================================================================

run_and_save_comparison <- function(comparison_name, file_prefix, sheet_name,
                                    group_vec, mat, batch_fac, purity_vec, sa_df,
                                    ds_out, ds_summary) {
    cat("\n  ---", comparison_name, "---\n")
    deg <- run_limma_group_comparison(mat, group_vec, batch_fac, purity_vec, sa_df)
    prefix_clean <- gsub("[^A-Za-z0-9_]", "_", file_prefix)

    if (!is.null(deg)) {
        fwrite(deg, file.path(ds_out, paste0("DPS_", prefix_clean, ".csv")))
        wb <- createWorkbook()
        addWorksheet(wb, sheet_name)
        writeData(wb, 1, deg)
        saveWorkbook(wb, file.path(ds_out, paste0("DPS_", prefix_clean, ".xlsx")), overwrite = TRUE)
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

    # --- Load phosphosite matrix with protein normalization ---
    mat0 <- tryCatch(
        load_phosphosite_matrix(ds_dir, protein_adjust = TRUE),
        error = function(e) {
            cat("  ERROR loading phosphosite matrix:", conditionMessage(e), "\n")
            NULL
        }
    )
    if (is.null(mat0)) { cat("  Skipping dataset\n\n"); next }

    # Impute and filter
    mat <- impute_and_filter_phospho(mat0, min_frac = min_frac_complete)
    cat("  After filtering:", nrow(mat), "phosphosites x", ncol(mat), "samples\n")

    # --- Get covariates ---
    bi <- get_batch_factor_phospho(ds_dir, colnames(mat))
    batch_fac <- if (!is.null(bi)) droplevels(bi$fac[colnames(mat)]) else NULL
    if (!is.null(bi)) cat("  Batch variable:", bi$name, "with", nlevels(batch_fac), "levels\n")

    purity_vec <- get_purity_covariate(ds_folder, ds_dir, colnames(mat))
    sa_df <- get_sex_age_covariates(ds_dir, colnames(mat))
    cat("  Covariates loaded (batch, purity, sex, age)\n")

    # --- Get TP53 classification ---
    tp53_ds <- tp53_all %>% filter(cancer_type == ds_cancer)
    tp53_mut <- tp53_ds %>% filter(mt == 1)

    wt_samples <- tp53_ds$sample_id[tp53_ds$wt == 1]
    gof_samples <- tp53_mut$sample_id[tp53_mut$GOF == 1]
    lof_samples <- tp53_mut$sample_id[tp53_mut$LOF == 1]
    hotspot_samples <- tp53_mut$sample_id[tp53_mut$hotspot == 1]
    dn_samples <- tp53_mut$sample_id[tp53_mut$DN == 1]
    nondn_samples <- tp53_mut$sample_id[tp53_mut$non_DN == 1]

    ds_out <- file.path(output_dir, ds_folder)
    dir.create(ds_out, recursive = TRUE, showWarnings = FALSE)
    ds_summary <- list(dataset = ds_folder, cancer_type = ds_cancer)

    # === 9 Comparisons ===

    # 1. TP53mt vs TP53wt
    g1 <- setNames(ifelse(tp53_ds$mt == 1, "TP53mt", "TP53wt"), tp53_ds$sample_id)
    g1 <- factor(g1[!is.na(g1)], levels = c("TP53wt", "TP53mt"))
    cat("    TP53wt:", sum(g1 == "TP53wt"), " TP53mt:", sum(g1 == "TP53mt"), "\n")
    ds_summary <- run_and_save_comparison("Comparison 1: TP53mt vs TP53wt", "TP53mt_vs_TP53wt", "mt_vs_wt",
        g1, mat, batch_fac, purity_vec, sa_df, ds_out, ds_summary)

    # 2. MUT_GOF vs MUT_LOF
    g2 <- setNames(c(rep("MUT_GOF", length(gof_samples)), rep("MUT_LOF", length(lof_samples))), c(gof_samples, lof_samples))
    g2 <- factor(g2, levels = c("MUT_LOF", "MUT_GOF"))
    ds_summary <- run_and_save_comparison("Comparison 2: GOF vs LOF", "MUT_GOF_vs_MUT_LOF", "GOF_vs_LOF",
        g2, mat, batch_fac, purity_vec, sa_df, ds_out, ds_summary)

    # 3. Hotspot vs MUT_LOF
    g3 <- setNames(c(rep("Hotspot", length(hotspot_samples)), rep("MUT_LOF", length(lof_samples))), c(hotspot_samples, lof_samples))
    g3 <- g3[!duplicated(names(g3))]; g3 <- factor(g3, levels = c("MUT_LOF", "Hotspot"))
    ds_summary <- run_and_save_comparison("Comparison 3: Hotspot vs LOF", "Hotspot_vs_MUT_LOF", "Hot_vs_LOF",
        g3, mat, batch_fac, purity_vec, sa_df, ds_out, ds_summary)

    # 4. MUT_GOF vs TP53wt
    g4 <- setNames(c(rep("MUT_GOF", length(gof_samples)), rep("TP53wt", length(wt_samples))), c(gof_samples, wt_samples))
    g4 <- factor(g4, levels = c("TP53wt", "MUT_GOF"))
    ds_summary <- run_and_save_comparison("Comparison 4: GOF vs TP53wt", "MUT_GOF_vs_TP53wt", "GOF_vs_wt",
        g4, mat, batch_fac, purity_vec, sa_df, ds_out, ds_summary)

    # 5. MUT_LOF vs TP53wt
    g5 <- setNames(c(rep("MUT_LOF", length(lof_samples)), rep("TP53wt", length(wt_samples))), c(lof_samples, wt_samples))
    g5 <- factor(g5, levels = c("TP53wt", "MUT_LOF"))
    ds_summary <- run_and_save_comparison("Comparison 5: LOF vs TP53wt", "MUT_LOF_vs_TP53wt", "LOF_vs_wt",
        g5, mat, batch_fac, purity_vec, sa_df, ds_out, ds_summary)

    # 6. Hotspot vs TP53wt
    g6 <- setNames(c(rep("Hotspot", length(hotspot_samples)), rep("TP53wt", length(wt_samples))), c(hotspot_samples, wt_samples))
    g6 <- factor(g6, levels = c("TP53wt", "Hotspot"))
    ds_summary <- run_and_save_comparison("Comparison 6: Hotspot vs TP53wt", "Hotspot_vs_TP53wt", "Hot_vs_wt",
        g6, mat, batch_fac, purity_vec, sa_df, ds_out, ds_summary)

    # 7. DN vs TP53wt
    g7 <- setNames(c(rep("DN", length(dn_samples)), rep("TP53wt", length(wt_samples))), c(dn_samples, wt_samples))
    g7 <- factor(g7, levels = c("TP53wt", "DN"))
    ds_summary <- run_and_save_comparison("Comparison 7: DN vs TP53wt", "DN_vs_TP53wt", "DN_vs_wt",
        g7, mat, batch_fac, purity_vec, sa_df, ds_out, ds_summary)

    # 8. Non-DN vs TP53wt
    g8 <- setNames(c(rep("NonDN", length(nondn_samples)), rep("TP53wt", length(wt_samples))), c(nondn_samples, wt_samples))
    g8 <- factor(g8, levels = c("TP53wt", "NonDN"))
    ds_summary <- run_and_save_comparison("Comparison 8: NonDN vs TP53wt", "NonDN_vs_TP53wt", "NonDN_vs_wt",
        g8, mat, batch_fac, purity_vec, sa_df, ds_out, ds_summary)

    # 9. DN vs non-DN
    g9 <- setNames(c(rep("DN", length(dn_samples)), rep("NonDN", length(nondn_samples))), c(dn_samples, nondn_samples))
    g9 <- g9[!duplicated(names(g9))]; g9 <- factor(g9, levels = c("NonDN", "DN"))
    ds_summary <- run_and_save_comparison("Comparison 9: DN vs NonDN", "DN_vs_NonDN", "DN_vs_NonDN",
        g9, mat, batch_fac, purity_vec, sa_df, ds_out, ds_summary)

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
fwrite(summary_df, file.path(output_dir, "DPS_summary_statistics.csv"))

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
saveWorkbook(wb_sum, file.path(output_dir, "DPS_summary_statistics.xlsx"), overwrite = TRUE)

cat("\n=== Differential Phosphosite Summary (padj < 0.05) ===\n\n")
print(as.data.frame(summary_df), row.names = FALSE)

cat("\n====================================================================\n")
cat("Pipeline Complete!\n")
cat("All results saved to:", output_dir, "\n")
cat("====================================================================\n")
