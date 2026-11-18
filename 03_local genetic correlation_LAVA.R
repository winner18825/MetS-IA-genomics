# =============================================================================
# 03 | Local genetic correlation -- LAVA
# =============================================================================
# Local Analysis of [co]Variant Association (LAVA: Werme J et al.,Nat Genet 2022)
# Runs LAVA (Local Analysis of [co]Variant Association) for each MetS-related
# trait against IA. For every pair the script:
#   1. standardises the GWAS summary statistics,
#   2. processes the 1000 Genomes EUR reference genotype data,
#   3. tests local univariate heritability and bivariate genetic correlation
#      locus by locus,
#   4. applies Benjamini-Hochberg FDR correction within each trait pair.
# =============================================================================

source("config.R")
library(data.table)
library(LAVA)
library(parallel)

# Paths and trait information

G1000_DIR <- file.path(REF_DIR, "g1000")
out_dir <- file.path(RESULTS_DIR, "LAVA")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ALL_TRAITS <- c(MT_TRAITS, IA_TRAIT)

# GWAS file paths
GWAS_FILES <- c(
  "FBG"          = file.path(DATA_DIR, "FBG.txt"),
  "HDL-C"        = file.path(DATA_DIR, "HDL-C.txt"),
  "Hypertension" = file.path(DATA_DIR, "Hypertension.txt"),
  "TG"           = file.path(DATA_DIR, "TG.txt"),
  "WC"           = file.path(DATA_DIR, "WC.txt"),
  "IA"           = file.path(DATA_DIR, "IA.txt")
)

# Case and control counts
# Continuous traits are set to NA.
trait_info <- data.table(phenotype = ALL_TRAITS, cases = NA_real_,
                         controls = NA_real_)
trait_info[phenotype == "Hypertension",`:=`(cases = 54358, controls = 408652)]
trait_info[phenotype == "IA",`:=`(cases = 6252, controls = 59544)]

# GWAS data standardisation
# Original columns: SNP, chr, pos, effect_allele, other_allele, 
#                   eaf, beta, se, pval, N, Z
# Output columns: SNP, A1, A2, N, Z
prepare_lava_sumstats <- function(infile, outfile) {
  df <- fread(infile)
  # Clean column names.
  setnames(df, trimws(names(df)))
  required <- c("SNP", "effect_allele", "other_allele", "N", "Z")
  missing_cols <- setdiff(required, names(df))
  df[, SNP := trimws(as.character(SNP))]
  df[, A1 := toupper(trimws(as.character(effect_allele)))]
  df[, A2 := toupper(trimws(as.character(other_allele)))]
  df[, N := suppressWarnings(as.numeric(N))]
  df[, Z := suppressWarnings(as.numeric(Z))]
  # Keep only the five columns LAVA requires.
  out <- df[, .(SNP, A1, A2, N, Z)]
  # Remove invalid records.
  out <- out[
    !is.na(SNP) & SNP != "" &
    !is.na(A1) & !is.na(A2) &
    A1 != A2 & !is.na(N) &
    N > 0 & !is.na(Z)
  ]
  out <- out[!duplicated(SNP)]
  fwrite(out, outfile, sep = "\t", na = "NA")
  invisible(out)
}

# Batch-generate six standardised GWAS files
GWAS_OUT_FILES <- setNames(
  file.path(DATA_DIR,paste0(ALL_TRAITS, "_lava.txt")),
  ALL_TRAITS
)

for (trait in ALL_TRAITS) {
  prepare_lava_sumstats(
    infile = GWAS_FILES[[trait]],
    outfile = GWAS_OUT_FILES[[trait]]
  )
}

# Locate the 1000 Genomes EUR reference genotype files
bed_files <- list.files(
  G1000_DIR, pattern = "\\.bed$",
  recursive = TRUE, full.names = TRUE, ignore.case = TRUE
)

idx <- which(
  tolower(basename(bed_files)) == "g1000_eur.bed"
)

if (length(idx) == 1) {
  bed_file <- bed_files[idx]
} else if (length(bed_files) == 1) {
  bed_file <- bed_files[1]
} else {
  stop(
    "Found multiple BED files or could not locate the G1000 EUR file; ",
    "please check the reference genotype files."
  )
}

# Reference file prefix.
REF_PREFIX <- sub("\\.bed$","",bed_file,ignore.case = TRUE)
required_ref <- paste0(REF_PREFIX,c(".bed", ".bim", ".fam"))

# Read the LAVA genome block file (GRCh37/hg19)
LOC_FILE <- file.path(G1000_DIR, "blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
loci <- read.loci(LOC_FILE)

SAMPLE_OVERLAP_FILE <- NULL

# Run the five LAVA analyses in parallel
N_CORES <- 5L  # 3 cores suggested by default; use 5 when memory allows.

run_one_pair <- function(MT) {
  IA <- "IA"
  cat("\nStarting analysis: ", MT, " vs ", IA, "\n")

  # Generate the input.info file.
  info <- data.table::copy(trait_info[match(c(MT, IA), phenotype)])
  info[, filename := basename(unname(GWAS_OUT_FILES[phenotype]))]
  data.table::setcolorder(info, c("phenotype", "cases", "controls", "filename"))

  INFO_FILE <- file.path(PROJECT_DIR, paste0("input.", MT, "_", IA, ".info.txt"))
  data.table::fwrite(info, INFO_FILE, sep = "\t", na = "NA")

  # Read the GWAS files and align SNPs.
  input <- tryCatch(
    process.input(input.info.file = INFO_FILE,
                  sample.overlap.file = SAMPLE_OVERLAP_FILE,
                  ref.prefix = REF_PREFIX,
                  phenos = c(MT, IA),
                  input.dir = DATA_DIR),
    error = function(e) {
      message(MT, " process.input failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(input)) {
    return(list(trait = MT, status = "input_failed", univ_n = 0, bivar_n = 0))
  }

  pair_univ_results <- list()
  pair_bivar_results <- list()
  failed_loci <- 0L

  # Iterate over all genome blocks.
  for (i in seq_len(nrow(loci))) {
    if (i %% 100 == 0) {
      message(MT, ": block ", i, "/", nrow(loci))
    }
    locus <- tryCatch(
      process.locus(loci[i, , drop = FALSE], input),
      error = function(e) {
        failed_loci <<- failed_loci + 1L
        NULL
      }
    )
    if (is.null(locus)) next
    res <- tryCatch(
      run.univ.bivar(locus, univ.thresh = 0.05),
      error = function(e) {
        failed_loci <<- failed_loci + 1L
        NULL
      }
    )
    if (is.null(res)) next
    # Block-level information.
    loc_info <- data.frame(
      MT_TRAIT = MT, IA_TRAIT = IA,
      LOC = loci$LOC[i], CHR = loci$CHR[i],
      START = loci$START[i], STOP = loci$STOP[i],
      n_snps = locus$n.snps, n_pcs = locus$K
    )

    # Save local heritability.
    if (!is.null(res$univ) && nrow(res$univ) > 0) {
      u <- as.data.frame(res$univ)
      pair_univ_results[[length(pair_univ_results) + 1L]] <- cbind(
        loc_info[rep(1L, nrow(u)), , drop = FALSE], u
      )
    }

    # Save local genetic correlation.
    if (!is.null(res$bivar) && nrow(res$bivar) > 0) {
      b <- as.data.frame(res$bivar)
      pair_bivar_results[[length(pair_bivar_results) + 1L]] <- cbind(
        loc_info[rep(1L, nrow(b)), , drop = FALSE], b
      )
    }
  }
  
  # Save local heritability results
  univ_n <- 0L
  bivar_n <- 0L
  significant_n <- NA_integer_

  if (length(pair_univ_results) > 0) {
    univ_all <- data.table::rbindlist(pair_univ_results, fill = TRUE)
    univ_n <- nrow(univ_all)

    data.table::fwrite(
      univ_all,
      file.path(out_dir, paste0(MT,"_", IA, "_LAVA_univariate_results.tsv")),
      sep = "\t"
    )
  }

  # Save local genetic correlation results
  # BH-FDR correction is applied within each trait pair.
  if (length(pair_bivar_results) > 0) {
    bivar_all <- data.table::rbindlist(pair_bivar_results, fill = TRUE)
    bivar_n <- nrow(bivar_all)

    if ("p" %in% names(bivar_all)) {
      bivar_all[, FDR_BH_pair := p.adjust(p, method = "BH")]
      data.table::setorder(bivar_all, p)
      significant_n <- sum(bivar_all$FDR_BH_pair < 0.05, na.rm = TRUE)
    }

    data.table::fwrite(
      bivar_all,
      file.path(out_dir, paste0(MT,"_", IA, "_LAVA_bivariate_results.tsv")),
      sep = "\t"
    )
  }
  return(list(
    trait = MT,status = "completed",
    univ_n = univ_n,bivar_n = bivar_n,
    significant_n = significant_n,
    failed_loci = failed_loci
  ))
}

# Create the parallel cluster
cl <- parallel::makeCluster(N_CORES, type = "PSOCK")

# Load the required packages on each worker.
parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(data.table)
    library(LAVA)
  })
  data.table::setDTthreads(1)
  NULL
})

# Pass variables and functions to the workers.
parallel::clusterExport(
  cl,
  varlist = c(
    "run_one_pair", "MT_TRAITS", "trait_info",
    "GWAS_OUT_FILES", "DATA_DIR","out_dir", 
    "G1000_DIR", "SAMPLE_OVERLAP_FILE","loci"
  ),
  envir = .GlobalEnv
)

# Run the five analyses in parallel
parallel_results <- tryCatch(
  parallel::parLapply(cl, MT_TRAITS, run_one_pair),
  finally = {
    parallel::stopCluster(cl)
    cat("\nParallel cluster stopped.\n")
  }
)

# Print run status
summary_results <- data.table::rbindlist(
  lapply(parallel_results, as.data.frame),
  fill = TRUE
)

cat("\n========== Five LAVA analyses complete ==========\n")
print(summary_results)
cat("\nResults directory: ", out_dir, "\n")
cat("\nBH-FDR correction was applied within each trait pair (no global correction).\n")
