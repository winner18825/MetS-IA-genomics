# Shared genetic architecture between metabolic syndrome-related traits and intracranial aneurysm

Analysis code for the study **"Shared genetic architecture between metabolic syndrome-related
traits and intracranial aneurysm: integrative genomic and transcriptomic evidence."**

This repository contains the scripts needed to reproduce the genome-wide
cross-trait analyses reported in the paper, from publicly available GWAS
summary statistics.

---

## Analysis pipeline

The study integrates five sequential modules. Each script corresponds to one
module and writes its output to `results/`.

| Step | Script | Method | Tool used |
|------|--------|--------|-----------|
| 1 | `01_genetic_correlation_LDSC.R` | Genome-wide genetic correlation | LDSC |
| 2 | `02_genetic_correlation_HDL.R`  | Genome-wide genetic correlation | HDL |
| 3 | `03_local genetic correlation_LAVA.R`  | Local genetic correlation | LAVA |
| 4 | `04_cross_trait_CPASSOC.R`      | Cross-trait meta-analysis (SHet) | CPASSOC |
| 5 | `05_cross_trait_PLACO.R`        | Cross-trait meta-analysis | PLACO |
| 6 | `06_merge_CPASSOC_PLACO.R`      | Intersection of pleiotropic loci | — |
| 7 | `07_colocalization.R`          | Bayesian colocalization | coloc |
| 8 | `08_Gene-based association_MAGMA.R`    | Gene-based association | MAGMA |
| 9 | `09_SMR_HEIDI.R`               | Gene-level pleiotropy | SMR + HEIDI |
| 10 | `10_bidirectional_MR.R`       | Bidirectional two-sample MR | TwoSampleMR + MR-PRESSO |

---

## Reproducibility note

All analyses were performed using publicly available R packages and open-source
implementations of the corresponding methods. The scripts provided in this
repository are fully reproducible and rely only on openly accessible,
citable software tools.

Users should ensure that the required package versions, reference datasets,
and external software dependencies are installed before running the scripts.
The exact versions of R packages and computational environments used for the
analysis are provided to facilitate reproducibility.

For one method that the authors distribute as plain R source files rather
than as installable packages — **PLACO** — you must download
the official source and place it under `tools/` (see paths in the script
headers). The scripts follow the documented workflow of each method; please
verify the exact function argument names against the version you download,
and confirm that results are consistent before publication.

---

## Software requirements

**R (>= 4.1)** with packages:

```r
install.packages(c("data.table", "MASS", "Matrix", "compiler","coloc", 
                   "ggplot2", "dplyr", "LAVA","stringr", "future", 
                   "TwoSampleMR", "MRPRESSO", "ieugwasr", "parallel"))
remotes::install_github("zhenin/HDL/HDL")
```

**External command-line tools** (install once, then set the paths in
`config.R`):

| Tool | Purpose | Source |
|------|---------|--------|
| LDSC | Genetic correlation | https://github.com/bulik/ldsc |
| SMR  | Gene-level pleiotropy | https://yanglab.westlake.edu.cn/software/smr/ |
| MAGMA  | Local genetic correlation | https://ctg.cncr.nl/software/magma/ |
| PLINK 1.9 | LD operations | https://www.cog-genomics.org/plink/ |

**Author-distributed R source files** (download and place under `tools/`):

| Tool | Place at | Source |
|------|----------|--------|
| PLACO   | `tools/PLACO/PLACO.R`     | https://github.com/RayDebashree/PLACO |

---

## Input data

GWAS summary statistics are **not redistributed** in this repository. Their
sources and access details are listed in **Supplementary Table S1** of the
paper.

Place the harmonised files under `data/gwas/` using these naming conventions:

- `data/gwas/<TRAIT>.txt` — used by LDSC, CPASSOC, MAGMA, colocalization.
  Expected columns: `SNP, CHR, BP, effect_allele, other_allele, eaf, beta,
  se, pval, N` (rename via the `COL` list in `config.R` if yours differ).
- `data/gwas/<TRAIT>_MTAG.txt` — used by HDL and PLACO; must include a
  Z-score column.
- `data/gwas/ma/<TRAIT>.ma` — used by SMR (SMR `.ma` format).

Trait names (`<TRAIT>`) are defined in `config.R`.

LD reference panels and GTEx v8 eQTL BESD files go under `reference/`
(see `config.R` and the `09_SMR_HEIDI.R` header for the expected layout).

---

## How to run

1. Clone the repository and open it as your R working directory.
2. Edit `config.R` — set the data, reference and tool paths for your machine.
3. Download the external tools and the PLACO source files.
4. Run the scripts in numerical order:

   ```r
   source("01_genetic_correlation_LDSC.R")
   source("02_genetic_correlation_HDL.R")
   source("03_local genetic correlation_LAVA.R")
   source("04_cross_trait_CPASSOC.R")
   source("05_cross_trait_PLACO.R")
   source("06_merge_CPASSOC_PLACO.R")
   source("07_colocalization.R")
   source("08_Gene-based association_MAGMA.R")
   source("09_SMR_HEIDI.R")
   source("10_bidirectional_MR.R")
   ```

All output is written to `results/`.

---

## Method references

- **LDSC** — Bulik-Sullivan et al. *Nat Genet* 2015.
- **HDL** — Ning, Pawitan & Shen. *Nat Genet* 2020.
- **LAVA** — Werme J et al. *Nat Genet* 2022.
- **CPASSOC** — Li et al. *Methods Mol Biol* 2017.
- **PLACO** — Ray & Chatterjee. *PLoS Genet* 2020.
- **coloc** — Giambartolomei et al. *PLoS Genet* 2014.
- **MAGMA** — de Leeuw CA et al. *PLoS Comput Biol* 2015.
- **SMR / HEIDI** — Zhu et al. *Nat Genet* 2016.
- **TwoSampleMR** — Hemani et al. *eLife* 2018.
- **MR-PRESSO** — Verbanck et al. *Nat Genet* 2018.

## Citation

If you use this code, please cite the associated paper (citation to be added
upon publication).

## License

Released under the MIT License — see [`LICENSE`](LICENSE).
