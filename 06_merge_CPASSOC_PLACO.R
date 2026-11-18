# =============================================================================
# 06 | Merge CPASSOC + PLACO significant loci
# =============================================================================
# A locus is retained as pleiotropic only if it is genome-wide significant in
# BOTH CPASSOC and PLACO. Significant SNPs from the two methods are
# intersected for every MetS-trait x IA pair, and all pairs are then stacked
# into a single table that feeds the colocalization step.
# =============================================================================

source("config.R")
library(data.table)

cpassoc_dir <- file.path(RESULTS_DIR, "CPASSOC")
placo_dir   <- file.path(RESULTS_DIR, "placo")
out_dir     <- file.path(RESULTS_DIR, "merge_cpassoc_placo")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

merged_all <- list()

get_raw_p <- function(path, new_col) {
  dat <- data.table::fread(path)
  # Identify the P-value column automatically.
  normalized_names <- tolower(gsub("[^[:alnum:]]", "", names(dat)))
  p_idx <- which(normalized_names %in% c("p", "pval", "pvalue"))
  if (length(p_idx) != 1L) {
    stop("Cannot uniquely identify the P-value column in: ", path, ". Available columns: ", paste(names(dat), collapse = ", "))
  }
  p_col <- names(dat)[p_idx]
  # Use the same SNP cleaning rules as CPASSOC.
  dat <- dat[grepl("^rs", SNP)]
  dat[, SNP := sub(",.*$", "", SNP)]
  dat <- dat[!duplicated(SNP)]
  # Keep only SNP and its original P-value.
  out <- dat[, .(SNP, p_tmp = as.numeric(get(p_col)))]
  data.table::setnames(out, "p_tmp", new_col)
  return(out)
}

for (MT in MT_TRAITS) {
  for (IA in IA_TRAIT) {
    tryCatch({
      cp_file <- file.path(cpassoc_dir, paste0(MT, "_", IA, "_CPASSOC_SHet.csv"))
      if (!file.exists(cp_file)) next
      cp <- data.table::fread(cp_file)

      # Read the original GWAS files for the two traits.
      trait1_path <- file.path(DATA_DIR, paste0(MT, ".txt"))
      trait2_path <- file.path(DATA_DIR, paste0(IA, ".txt"))

      # Extract the original single-trait P-values.
      p1 <- get_raw_p(trait1_path, "pval_trait1")
      p2 <- get_raw_p(trait2_path, "pval_trait2")

      # Match the original P-values to the CPASSOC results by SNP.
      cp <- merge(cp, p1, by = "SNP")
      cp <- merge(cp, p2, by = "SNP")

      pl_file <- file.path(placo_dir, paste0(MT, "_", IA, "_PLACO_result.txt"))
      pl <- data.table::fread(pl_file)

      # CPASSOC: genome-wide significant SHet, both single-trait p < 1e-3.
      cp <- cp[p_SHet < 5e-8 & pval_trait1 < 1e-3 & pval_trait2 < 1e-3]
      # PLACO: genome-wide significant composite p.
      pl <- pl[p.placo < 5e-8]

      shared <- merge(cp, pl, by = "SNP", suffixes = c(".x", ".y"))
      if (nrow(shared) == 0) next
      shared[, prefix := paste0(MT, "&", IA)]
      data.table::fwrite(shared,
                         file.path(out_dir, paste0(MT, "_", IA, "_merge_cpassoc_placo_sig.txt")),
                         sep = "\t")
      merged_all[[paste0(MT, "_", IA)]] <- shared
    }, error = function(e) {
      message("ERROR: ", MT, " - ", IA, " : ", e$message)
      # Continue to the next pair.
    })
  }
}

all_loci <- data.table::rbindlist(merged_all, fill = TRUE)
data.table::fwrite(all_loci,
                   file.path(out_dir, "merged_all_placo_sig.txt"), sep = "\t")

cat("Merged", nrow(all_loci), "pleiotropic loci across",
    length(merged_all), "trait pairs.\n")
cat("Combined table:",
    file.path(out_dir, "merged_all_placo_sig.txt"), "\n")
