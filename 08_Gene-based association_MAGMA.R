# =============================================================================
# 08 | Gene-based association -- MAGMA
# =============================================================================
# Multi-marker Analysis of GenoMic Annotation (MAGMA: de Leeuw CA et al.,
# PLoS Comput Biol 2015)
# Runs MAGMA v1.10 gene-based analysis for the six traits, then intersects the
# gene-level results with the shared (CPASSOC + PLACO) lead-SNP loci to build
# candidate pleiotropic genes. A gene is retained only when it lies in the
# pair-specific +/- 500 kb shared region and is FDR-significant in both traits
# of the pair. BH-FDR correction is applied separately within each phenotype.
# =============================================================================

library(data.table)

GENE_LOC <- file.path(REF_DIR, "MAGMA/NCBI37.3.gene.loc")
LD_BFILE <- file.path(REF_DIR, "1000G_EUR/EUR")
LOCUS_FILE <- file.path(RESULTS_DIR, "merge_cpassoc_placo/merged_all_placo_sig.txt")
MAGMA_BIN <- file.path(MAGMA_ROOT, "magma")

out_dir <- file.path(RESULTS_DIR, "MAGMA")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

TRAITS <- c(MT_TRAITS, IA_TRAIT)
BINARY_TOTAL_N <- c(Hypertension = 54358 + 408652, IA = 6252 + 59544)
LOCUS_WINDOW <- 500000L
GENE_WINDOW_UP <- 10L
GENE_WINDOW_DOWN <- 10L

# Standardise chromosome labels.
norm_chr <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("^CHR", "", x)
  x <- sub("^0+([0-9]+)$", "\\1", x)
  x
}

# MAGMA runner
run_magma <- function(args, logfile) {
  output <- suppressWarnings(system2(MAGMA_BIN, args = args, stdout = TRUE, stderr = TRUE))
  writeLines(as.character(output), logfile)
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  if (status != 0L) stop("MAGMA run failed; check the log: ", logfile)
  invisible(output)
}

# Read the shared lead SNPs and their prefixes
loci <- fread(LOCUS_FILE)
setnames(loci, toupper(trimws(names(loci))))

loci[, SNP := trimws(as.character(SNP))]
loci[, CHR := norm_chr(CHR)]
loci[, BP := as.numeric(BP)]
loci[, PREFIX := gsub("\\s+", "", as.character(PREFIX))]

loci <- loci[
  !is.na(SNP) & SNP != "" &
    CHR %in% as.character(1:22) &
    is.finite(BP) & BP > 0 &
    !is.na(PREFIX) & PREFIX != ""
]

prefix_mtrts <- strsplit(loci$PREFIX, "&", fixed = TRUE)

# Ensure each group is one MetS trait mtired with IA.
mt_trait <- vapply(prefix_parts, function(x) {
  x <- trimws(x)
  if (sum(x == "IA") != 1) return(NA_character_)
  setdiff(x, "IA")
}, character(1))

loci[, TRAIT1 := mt_trait]
loci[, TRAIT2 := "IA"]
loci[, PREFIX := paste(TRAIT1, TRAIT2, sep = "&")]

loci <- unique(loci, by = c("PREFIX", "SNP", "CHR", "BP"))

loci[, locus_id := paste0(PREFIX, "_L", .I)]
loci[, LOC_START := pmax(1, BP - LOCUS_WINDOW)]
loci[, LOC_STOP := BP + LOCUS_WINDOW]

# Read the GRCh37 gene location file
gene_raw <- read.table(GENE_LOC, header = FALSE, comment.char = "#", 
                       stringsAsFactors = FALSE, fill = TRUE)

genes <- data.table(
  GENE_ID = as.character(gene_raw[[1]]),
  G_CHR = norm_chr(gene_raw[[2]]),
  G_START = as.numeric(gene_raw[[3]]),
  G_STOP = as.numeric(gene_raw[[4]])
)

genes <- genes[
  !is.na(GENE_ID) & GENE_ID != "" &
    G_CHR %in% as.character(1:22) &
    is.finite(G_START) & is.finite(G_STOP)
]

genes[G_START > G_STOP, c("G_START", "G_STOP") := .(G_STOP, G_START)]
gene_chr_check <- genes[, .(n_chr = uniqueN(G_CHR)), by = GENE_ID]

# Collapse multiple location records for the same gene ID.
genes <- genes[, .(
  G_CHR = G_CHR[1],
  G_START = min(G_START),
  G_STOP = max(G_STOP)
), by = GENE_ID]

# Extend each gene by 10 kb on each side.
genes[, G_START_EXT := pmax(1, G_START - GENE_WINDOW_UP * 1000)]
genes[, G_STOP_EXT := G_STOP + GENE_WINDOW_DOWN * 1000]

# Convert gene IDs to gene symbols.
genes[, GENE_SYMBOL := NA_character_]

if (requireNamespace("AnnotationDbi", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  tryCatch({
    db <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
    gene_symbols <- AnnotationDbi::mapIds(
      db, keys = genes$GENE_ID, column = "SYMBOL",
      keytype = "ENTREZID", multiVals = "first"
    )
    genes[, GENE_SYMBOL := unname(gene_symbols[GENE_ID])]
  }, error = function(e) message("Gene symbol conversion failed; keeping gene IDs: ", conditionMessage(e)))
}

# Build the shared-locus to candidate-gene mapping by prefix
lead_map <- loci[, .(
  PREFIX, TRAIT1, TRAIT2, locus_id, SNP,
  LEAD_CHR = CHR, LEAD_BP = BP,
  L_START = LOC_START, L_STOP = LOC_STOP
)]

genes_overlap <- genes[, .(
  GENE_ID, GENE_SYMBOL, G_START, G_STOP,
  G_CHR, G_START_EXT, G_STOP_EXT
)]

setkeyv(lead_map, c("LEAD_CHR", "L_START", "L_STOP"))
setkeyv(genes_overlap, c("G_CHR", "G_START_EXT", "G_STOP_EXT"))

locus_gene_map <- foverlaps(
  x = lead_map,
  y = genes_overlap,
  by.x = c("LEAD_CHR", "L_START", "L_STOP"),
  by.y = c("G_CHR", "G_START_EXT", "G_STOP_EXT"),
  type = "any",
  nomatch = 0L
)

# Distance from the lead SNP to the gene body.
locus_gene_map[, distance_bp := pmax(G_START - LEAD_BP, 0, LEAD_BP - G_STOP)]

locus_gene_map[, mapping_type := ifelse(
  distance_bp == 0, "gene_body",
  ifelse(distance_bp <= GENE_WINDOW_UP * 1000, "gene_flank_10kb", "within_500kb_locus")
)]

# Save the full candidate-gene mapping.
fwrite(locus_gene_map, 
       file.path(out_dir, "MAGMA_shared_locus_candidate_genes.tsv"), sep = "\t")

# Save candidate genes per prefix.
PREFIX_DIR <- file.path(out_dir, "by_prefix")
dir.create(PREFIX_DIR, recursive = TRUE, showWarnings = FALSE)

for (prefix in unique(loci$PREFIX)) {
  safe_prefix <- gsub("[^A-Za-z0-9_-]", "_", prefix)
  fwrite(loci[PREFIX == prefix], file.path(PREFIX_DIR, paste0(safe_prefix, "_lead_SNPs.tsv")), sep = "\t")
  fwrite(locus_gene_map[PREFIX == prefix], file.path(PREFIX_DIR, paste0(safe_prefix, "_candidate_genes.tsv")), sep = "\t")
}

# Read the six full GWAS and prepare the MAGMA input
snp_loc_list <- list()
for (trait in TRAITS) {
  infile <- file.path(DATA_DIR, paste0(trait, ".txt"))

  d <- fread(infile)
  setnames(d, toupper(trimws(names(d))))

  d[, SNP := trimws(as.character(SNP))]
  d[, CHR := norm_chr(CHR)]
  d[, BP := as.numeric(POS)]
  d[, P := as.numeric(PVAL)]
  d[, N := as.numeric(N)]

  d <- d[
    !is.na(SNP) & SNP != "" &
      CHR %in% as.character(1:22) &
      is.finite(BP) & BP > 0 &
      is.finite(P) & P >= 0 & P <= 1 &
      is.finite(N) & N > 0
  ]

  # Remove duplicate SNP IDs.
  dup <- d[, .N, by = SNP][N > 1, SNP]
  d <- d[!SNP %chin% dup]

  # P-values cannot be exactly zero.
  d[, P := pmax(P, .Machine$double.xmin)]

  # For binary traits, check that N is the total sample size.
  if (trait %in% names(BINARY_TOTAL_N)) {
    expected_n <- BINARY_TOTAL_N[[trait]]
    median_n <- median(d$N, na.rm = TRUE)
    if (median_n < expected_n * 0.8) {
      warning(trait, " median N is markedly below cases + controls; verify whether N is the total sample size.")
    }
  }

  # Generate the MAGMA p-value file.
  PVAL_FILE <- file.path(out_dir, paste0(trait, ".pval"))
  fwrite(d[, .(SNP, P, N)], PVAL_FILE, sep = "\t")

  # Save SNP positions.
  snp_loc_list[[trait]] <- d[, .(SNP, CHR, BP, GWAS_trait = trait)]
}

# Use the 1000G EUR BIM file as the SNP location file
SNP_LOC_FILE <- paste0(LD_BFILE, ".bim")
bim <- fread(SNP_LOC_FILE, header = FALSE)
setnames(bim, c("CHR", "SNP", "CM", "BP", "A1", "A2"))

# Check SNP ID overlap between each GWAS and the reference.
ref_snps <- unique(bim$SNP)

for (trait in TRAITS) {
  pval_file <- file.path(out_dir, paste0(trait, ".pval"))
  d <- fread(pval_file)
  matched <- sum(d$SNP %chin% ref_snps)
  cat(trait, " matched SNPs: ", matched, "/", nrow(d),
      " (", round(matched / nrow(d) * 100, 2), "%)\n")
}

# Run the MAGMA SNP-to-gene annotation once
ANNOT_PREFIX <- file.path(out_dir, "all_traits_annotation")

run_magma(
  c("--annotate",
    paste0("window=", GENE_WINDOW_UP, ",", GENE_WINDOW_DOWN),
    "--snp-loc", shQuote(SNP_LOC_FILE),
    "--gene-loc", shQuote(GENE_LOC),
    "--out", shQuote(ANNOT_PREFIX)),
  file.path(out_dir, "annotation.log")
)

GENE_ANNOT <- paste0(ANNOT_PREFIX, ".genes.annot")

# Run MAGMA gene-based analysis for the six traits
gene_results <- list()

for (trait in TRAITS) {
  PVAL_FILE <- file.path(out_dir, paste0(trait, ".pval"))
  OUT_PREFIX <- file.path(out_dir, paste0(trait, "_MAGMA"))

  tryCatch({
    run_magma(
      c(
        "--bfile", shQuote(LD_BFILE),
        "--pval", shQuote(PVAL_FILE),
        "use=SNP,P",
        "ncol=N",
        "--gene-annot", shQuote(GENE_ANNOT),
        "--out", shQuote(OUT_PREFIX)
      ),
      file.path(out_dir, paste0(trait, "_MAGMA.log"))
    )
    gene_file <- paste0(OUT_PREFIX, ".genes.out")

    gr <- fread(gene_file)
    setnames(gr, toupper(trimws(names(gr))))

    gr[, GENE_ID := as.character(GENE)]
    gr[, P := as.numeric(P)]
    gr[, trait := trait]

    gene_results[[trait]] <- gr

  }, error = function(e) {
    message(trait, " MAGMA analysis failed: ", conditionMessage(e))
  })
}

# Combine the six MAGMA results and apply BH-FDR per trait

all_gene_results <- rbindlist(gene_results, fill = TRUE)
gene_symbol_map <- unique(genes[, .(GENE_ID, GENE_SYMBOL)])

all_gene_results <- merge(
  all_gene_results,
  gene_symbol_map,
  by = "GENE_ID",
  all.x = TRUE,
  sort = FALSE
)

# BH-FDR correction separately within each of the six traits.
all_gene_results[, FDR_BH := p.adjust(P, method = "BH"), by = trait]

fwrite(all_gene_results,
       file.path(out_dir, "MAGMA_all_gene_results.tsv"),sep = "\t")

# Extract, per prefix, the MAGMA results of the two traits
pair_info <- unique(loci[, .(PREFIX, TRAIT1, TRAIT2)])

pair_gene_results <- rbindlist(lapply(seq_len(nrow(pair_info)), function(i) {

  prefix <- pair_info$PREFIX[i]
  t1 <- pair_info$TRAIT1[i]
  t2 <- pair_info$TRAIT2[i]

  # Candidate genes for the current prefix.
  candidate <- unique(locus_gene_map[PREFIX == prefix, .(
    PREFIX, TRAIT1, TRAIT2, locus_id, SNP,
    LEAD_CHR, LEAD_BP, GENE_ID, distance_bp, mapping_type
  )])

  if (nrow(candidate) == 0) return(NULL)

  # Restrict to the two traits of the current prefix.
  gr <- all_gene_results[trait %chin% c(t1, t2)]

  if (nrow(gr) == 0) return(NULL)

  merge(candidate, gr, by = "GENE_ID",
        all = FALSE,allow.cartesian = TRUE)
}), fill = TRUE)

fwrite(pair_gene_results,
       file.path(out_dir, "MAGMA_pairwise_gene_results.tsv"),sep = "\t")

# Filter MAGMA-significant genes per prefix
# Remove duplicate records of the same prefix/gene/trait across lead SNPs.
pair_gene_results_unique <- unique(
  pair_gene_results,by = c("PREFIX", "GENE_ID", "trait")
)

pair_gene_sig <- pair_gene_results_unique[
  !is.na(FDR_BH) & FDR_BH < 0.05
]

fwrite(pair_gene_sig,
       file.path(out_dir, "MAGMA_pairwise_significant_genes.tsv"),sep = "\t")

# Keep genes significant in both traits of a prefix
pair_shared <- pair_gene_sig[, .(
  n_traits = uniqueN(trait),
  significant_traits = paste(sort(unique(trait)), collapse = "&")
), by = .(PREFIX, GENE_ID)]

# Both traits must be significant.
pair_shared <- pair_shared[n_traits == 2]

pair_shared_results <- merge(
  pair_gene_sig, pair_shared,
  by = c("PREFIX", "GENE_ID"), all = FALSE)

fwrite(pair_shared_results,
       file.path(out_dir, "MAGMA_pairwise_both_traits_significant.tsv"),sep = "\t")

# Save per-prefix results
for (prefix in pair_info$PREFIX) {
  safe_prefix <- gsub("[^A-Za-z0-9_-]", "_", prefix)

  tmp_all <- pair_gene_results[PREFIX == prefix]
  tmp_sig <- pair_gene_sig[PREFIX == prefix]
  tmp_both <- pair_shared_results[PREFIX == prefix]

  if (nrow(tmp_all) > 0) {
    fwrite(tmp_all, file.path(PREFIX_DIR, paste0(safe_prefix, "_MAGMA_all.tsv")), sep = "\t")
  }
  if (nrow(tmp_sig) > 0) {
    fwrite(tmp_sig, file.path(PREFIX_DIR, paste0(safe_prefix, "_MAGMA_significant.tsv")), sep = "\t")
  }
  if (nrow(tmp_both) > 0) {
    fwrite(tmp_both, file.path(PREFIX_DIR, paste0(safe_prefix, "_MAGMA_both_traits_significant.tsv")), sep = "\t")
  }
}

# Summarise per-prefix results
prefix_summary <- pair_gene_results_unique[, .(
  Candidate_genes = uniqueN(GENE_ID),
  Significant_gene_trait = sum(FDR_BH < 0.05, na.rm = TRUE),
  Significant_genes = uniqueN(GENE_ID[FDR_BH < 0.05])
), by = PREFIX]

both_summary <- pair_shared[, .(
  Both_traits_significant_genes = .N
), by = PREFIX]

prefix_summary <- merge(prefix_summary, both_summary, by = "PREFIX", all.x = TRUE)
prefix_summary[is.na(Both_traits_significant_genes), Both_traits_significant_genes := 0L]

fwrite(prefix_summary,
       file.path(out_dir, "MAGMA_prefix_summary.tsv"),sep = "\t")
