# Population-stratified archaic introgression and selection

This folder contains the complete analysis code for investigating **population-stratified archaic introgression and selection signals** in the human genome. The code is organized following the order of results presented in the corresponding section of our study, covering data processing, statistical analysis, and visualization of introgression patterns, haplotype structures, geographic distribution, extended haplotype homozygosity (EHH), and expression associations.

## Folder Contents and Functionality

### 1. Gene content and regulatory annotation of introgressed archaic chunks  
**File:** `01.Gene_content_and_regulatory_annotation_of_introgressed_archaic_chunks.Rmd`  
- **Purpose:** Annotates introgressed archaic genomic blocks with gene features (exons, introns, UTRs) and regulatory elements (cCREs, TF binding sites) using Ensembl GFF3 and ENCODE cCRE data.  
- **Key outputs:**  
  - Gene track visualizations highlighting introgressed sites and regulatory variants.  
  - Excel/TSV files listing introgressed variants with their functional consequences.  
- **Dependencies:** `rtracklayer`, `ggplot2`, `ggrepel`, `dplyr`, `tidyr`.

### 2. Haplotype heatmap of the introgression region  
**File:** `02.haplotype_heatmap_of_the_introgression_region.Rmd`  
- **Purpose:** Generates heatmaps of haplotype patterns across introgressed regions, with sample annotations by ancestry.  
- **Steps:**  
  - Filters SNPs based on frequency and missingness.  
  - Prepares genotype matrix for `pheatmap`.  
  - Plots clustered heatmaps with optional column annotations (introgression sites, variant types).  
- **Key outputs:** PDF heatmaps with informative file names (region length, variant count).  
- **Dependencies:** `pheatmap`, `RColorBrewer`, `dplyr`.

### 3. Haplotype networks  
**File:** `03.Haplotype_networks.Rmd`  
- **Purpose:** Constructs haplotype networks from phased genotype data using the `geneHapR` package, and calculates pairwise differences between haplotypes.  
- **Workflow:**  
  - `table2hap` → `hap_summary` → `get_hapNet`.  
  - Pairwise difference calculation and visualization (heatmaps, networks).  
  - Optional filtering to retain only common haplotypes + archaic haplotype.  
- **Key outputs:** Haplotype summary tables, difference matrices, network plots, analysis reports.  
- **Dependencies:** `geneHapR`, `igraph`, `pheatmap`, `tidyverse`.

### 4. Geographic distribution of archaic alleles (1KGP + HGDP)  
**File:** `04.Geographic_distribution_of_ArchaicAllele_1KGP_HGDP_data.Rmd`  
- **Purpose:** Maps the frequency of introgressed archaic alleles across global populations using HGDP and 1KGP data.  
- **Key features:**  
  - Parses VCF files, converts genotypes, and merges with population metadata.  
  - Calculates allele frequencies per population and ancestry group.  
  - Plots world maps with pie charts showing archaic vs. non-archaic allele frequencies.  
  - Splits East/Southeast Asian populations if desired.  
- **Key outputs:** PDF maps, frequency tables by population and ancestry.  
- **Dependencies:** `scatterpie`, `ggrepel`, `cowplot`, `rworldmap` (via `ggplot2::borders`).

### 5. EHH analysis – file preparation  
**File:** `05.EHH_analysis_01_file_prepare.sh`  

- **Purpose:** Bash script to prepare input files for extended haplotype homozygosity (EHH) analysis.  
- **Steps:**  
  - Extracts a focal region around a marker using `bcftools`.  
  - Filters samples by population (EAS, CHB, CHS).  
  - Calls an R script (`../ehh.R`) to compute EHH.  
- **Input:** A `demo.info` file with chromosome, position, and marker ID.  
- **Output:** Per-population EHH PDF plots.

### 6. EHH analysis – plotting  
**File:** `05.EHH_analysis_02_plot.Rmd`  
- **Purpose:** R script (called by the bash script) to compute and plot EHH using the `rehh` package.  
- **Key function:** `data2haplohh` + `calc_ehh` + `plot`.  
- **Dependencies:** `rehh`, `argparse`.

### 7. PRDM16 allele frequency vs. latitude correlation  
**File:** `06.PRDM16_correlation_alleleFrequency_Latitude.Rmd`  
- **Purpose:** Tests correlation between introgressed allele frequency and latitude for functional sites (e.g., regulatory variants) in PRDM16.  
- **Workflow:**  
  - Loops over candidate sites, reads precomputed frequency files.  
  - Performs Pearson correlation with latitude (filtering for East Asian populations).  
  - Generates scatter plots with regression lines and correlation statistics.  
- **Key outputs:** PDF correlation plots, summary statistics table.  
- **Dependencies:** `ggpubr`, `scales`, `tidyverse`.

### 8. Expression association analysis – APG data  
**File:** `07.expression.APG.py`  
- **Purpose:** Tests association between archaic allele genotype and gene expression levels in the **APG** dataset.  
- **Input:**  
  - SV/SNP list (from `rg -f test.list Nean_SV.pvalue`)  
  - Expression matrix (RNA_APG_MGI.OUT.gene_counts.TMM)  
  - Gene ID (e.g., `ENSG00000169155.10`)  
  - VCF file with phased genotypes  
  - Output prefix  
- **Workflow:**  
  - Pre-loads expression data for the target gene.  
  - For each variant, extracts genotypes from VCF, classifies samples as `intro/intro`, `intro/non-intro`, `non-intro/non-intro`.  
  - Performs OLS regression of expression on genotype.  
  - Generates boxplots + swarmplots and saves to a multi-page PDF.  
- **Output:** PDF with one page per variant, showing expression distribution by genotype.  
- **Dependencies:** `pysam`, `pandas`, `matplotlib`, `seaborn`, `statsmodels`.

### 9. Expression association analysis – MAGE data  
**File:** `07.expression.MAGE.py`  
- **Purpose:** Same as above, but for the **MAGE** dataset (different expression matrix format and VCF path).  
- **Differences:**  
  - Expression file is `TMM.filtered.TSS.MAGE.v1.0.bed`.  
  - VCF path is hardcoded (`1KGP.CHM13v2.0...`).  
  - Requires only 4 arguments (no expression file path, as it's fixed).  
- **Usage example:**  
  
  ```bash
  python expression.MAGE.2.py <(rg -f test.list Nean_SV.pvalue) ENSG00000169155.10 ZBTB43.MAGE
  ```
- **Output:** Multi-page PDF with expression-by-genotype plots.

## How to Use the Code

### General Requirements
- **R version ≥ 4.0** with packages: `tidyverse`, `pheatmap`, `RColorBrewer`, `geneHapR`, `igraph`, `scatterpie`, `ggrepel`, `cowplot`, `rehh`, `argparse`, `rtracklayer`, `ggpubr`, `scales`.
- **Python 3** with: `pysam`, `pandas`, `matplotlib`, `seaborn`, `statsmodels`, `numpy`.
- **External tools:** `bcftools`, `tabix`, `rg` (ripgrep, optional).

### Running the R Notebooks
Each `.Rmd` file is self-contained and can be knitted in RStudio or run via `rmarkdown::render()`.  
**Important:** Before running, update the file paths at the top of each notebook to point to your local data directories (population info, VCF files, expression matrices, etc.). The scripts use absolute paths; you may need to modify them.

### Running the EHH Analysis (bash + R)
1. Prepare a `demo.info` file with three tab-separated columns:  
   ```
   chr   position   markerID
   ```
   Example:  
   ```
   chr1  3320163   PRDM16_snp1
   ```
2. Run the preparation script:  
   ```bash
   bash 05.EHH_analysis_01_file_prepare.sh
   ```
   This will generate population-specific EHH plots (`*_EAS.pdf`, `*_CHB.pdf`, `*_CHS.pdf`).

### Running the Expression Association Scripts
#### APG version
```bash
python 07.expression.APG.py <(rg -f test.list Nean_SV.pvalue) RNA_APG_MGI.OUT.gene_counts.TMM ENSG00000169155.10 ../05.Blocks/Nean-gene/Nean_4229/Nean_4229.All.vcf.gz ZBTB43.MGI
```
- `<(rg -f test.list Nean_SV.pvalue)` supplies a list of variants (chrom, pos, etc.) via process substitution.  
- The expression file is tab-separated with a header; first column is gene ID, subsequent columns are samples.  
- The VCF must be indexed (`bcftools index`).  
- Output PDF: `ZBTB43.MGI.pdf`.

#### MAGE version
```bash
python 07.expression.MAGE.2.py <(rg -f test.list Nean_SV.pvalue) ENSG00000169155.10 ZBTB43.MAGE
```
- The expression file (`TMM.filtered.TSS.MAGE.v1.0.bed`) and VCF path are hardcoded in the script; adjust if needed.  
- Output PDF: `ZBTB43.MAGE.pdf`.

## Notes
- All scripts assume a specific directory structure. Modify paths according to your setup.
- The R notebooks often use `conflicted` to resolve function name clashes; ensure it is installed.

## Contact
For questions or issues, please contact the corresponding author or open an issue on the GitHub repository.

