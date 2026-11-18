# =============================================================================
# config.R  --  central configuration
# =============================================================================
# Edit the paths and parameters below ONCE to match your environment.
# Every analysis script begins with  source("config.R")  and inherits them.
# Run all scripts with the repository root as the working directory.
# =============================================================================

## ---- Project directories ----------------------------------------------------
# Folder holding the harmonised GWAS summary statistics.
# Naming convention (kept identical to the original analysis):
#   <TRAIT>.txt        used by LDSC, CPASSOC, MAGMA and colocalization
#   <TRAIT>_MTAG.txt   used by HDL and PLACO (must contain a Z-score column)
#   ma/<TRAIT>.ma      used by SMR  (SMR .ma format: SNP A1 A2 freq b se p n)
DATA_DIR    <- "data"      # GWAS summary statistics (NOT redistributed)
REF_DIR     <- "reference"      # LD reference panels, eQTL BESD files, etc.
RESULTS_DIR <- "results"        # all analysis output is written here

## ---- Trait definitions ------------------------------------------------------
MT_TRAITS <- c("FBG","HDL-C","Hypertension","TG","WC")
IA_TRAIT <- "IA"
## ---- Expected GWAS column names ---------------------------------------------
# The <TRAIT>.txt files are expected to contain these columns. Adjust the
# values if your files use different headers.
COL <- list(
  SNP = "SNP", CHR = "chr", BP = "pos",
  EA  = "effect_allele", OA = "other_allele",
  EAF = "eaf", BETA = "beta", SE = "se", P = "pval",
  N="N"
)

## ---- External command-line tools --------------------------------------------
# LDSC  -- https://github.com/bulik/ldsc  (Python 2.7 command-line tool)
LDSC_DIR <- "tools/ldsc"        # folder containing ldsc.py / munge_sumstats.py
LDSC_PY  <- "python2"           # interpreter used to run LDSC
# 
# SMR   -- https://yanglab.westlake.edu.cn/software/smr/
SMR_BIN  <- "tools/smr/smr"     # path to the smr executable

# PLINK 1.9 -- https://www.cog-genomics.org/plink/
PLINK_BIN <- "tools/plink/plink"

# MAGMA 1.10 -- https://ctg.cncr.nl/software/magma/
MAGMA_ROOT <- "tools/MAGMA/MAGMA_bin"

## ---- Create the results directory if missing --------------------------------
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)