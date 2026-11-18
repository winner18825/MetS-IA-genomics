# =============================================================================
# 04 | Cross-trait meta-analysis -- CPASSOC (SHet)
# =============================================================================
# Detects pleiotropic loci shared by MetS-related traits and IA 
# using the SHet statistic from CPASSOC (Li et al., Methods Mol Biol 2017).
# The SHet statistic accounts for cross-trait association signals by jointly
# analyzing summary statistics and is designed to accommodate heterogeneity
# of effect sizes across phenotypes.
#
# The CPASSOC functions used in this analysis were implemented based on the
# original method description and workflow. The required functions are provided
# locally in: CPASSOC_FunctionSet.R
# =============================================================================

source("config.R")
source("CPASSOC_FunctionSet.R")
library(dplyr)
library(stringr)
library(data.table)
library(future)

out_dir <- file.path(RESULTS_DIR, "CPASSOC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_cpassoc <- function(MT, IA) {
  trait1_path <- file.path(DATA_DIR, paste0(MT, ".txt"))
  trait2_path <- file.path(DATA_DIR, paste0(IA, ".txt"))

  trait1_data <- fread(trait1_path, data.table = FALSE)
  trait2_data <- fread(trait2_path, data.table = FALSE)

  trait1_samplesize <- max(trait1_data$N, na.rm = TRUE)
  trait2_samplesize <- max(trait2_data$N, na.rm = TRUE)

  sample_size <- c(trait1_samplesize, trait2_samplesize)
  d <- cpassoc.combine_assoc(files_path = c(trait1_path, trait2_path), traits = c(MT,IA))

  # Restrict to independent SNPs and derive the cross-trait Z-score correlation.
  cor_mat <- d %>%
    dplyr::filter(SNP %in% snplist) %>%
    dplyr::select(tidyselect::starts_with("z")) %>%
    as.matrix() %>%
    cor()

  # CPASSOC SHet/SHom tests.
  sh <- cpassoc.shet_shom_test(
    gwas_data = d, SHet_test = TRUE, SHom_test = TRUE,
    workers = 4, mvrnorm_n = 1e6, set_seed = 666)

  # Save results.
  fwrite(sh, file = paste0(out_dir, "/", MT, "_", IA, "_CPASSOC_SHet.csv"),
         quote = FALSE, sep = "\t")
}

for (MT in MT_TRAITS) {
  gwas_dat <- data.table::fread(file.path(DATA_DIR, paste0(MT, ".txt")), data.table = FALSE)

  # Use PLINK to select a set of independent SNPs.
  fn <- tempfile()
  write.table(data.frame(SNP = gwas_dat[["SNP"]], P = gwas_dat[["pval"]]),
              file = fn, row.names = F, col.names = T, quote = F)
  plink_path <- file.path(PLINK_BIN, "plink_linux", "plink")
  args <- c("--bfile", file.path(REF_DIR, "1kg.v3/EUR"),
            "--clump", fn,
            "--clump-p1", 5e-8,
            "--clump-r2", 0.1,
            "--clump-kb", 1000,
            "--out", "plink_clump")
  system2(command = plink_path, args = args)

  # Remove the temporary file.
  unlink(fn)

  # Read plink_clump.clumped and obtain the set of independent SNPs.
  snplist_clump <- data.table::fread("plink_clump.clumped", data.table = F)
  snplist1 <- snplist_clump$SNP
  snplist2 <- unlist(strsplit(gsub("\\([^()]*\\)", "", snplist_clump$SP2), split = ","))
  snplist <- unique(c(snplist1, snplist2))

  for (IA in IA_TRAIT) {
    result <- run_cpassoc(MT, IA)
  }
}

cat("CPASSOC cross-trait meta-analysis complete. Results in:", out_dir, "\n")
