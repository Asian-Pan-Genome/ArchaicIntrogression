# Three Denisovan introgression components in East Eurasian genomes



# Gaussian Mixture Modeling (GMM) pipeline

This document describes the R-based workflow for inferring the optimal number of Denisovan introgression components per haploid genome using Gaussian mixture modeling (GMM) on match rates (MR) of introgressed segments. The pipeline implements likelihood‑ratio tests (LRT) to compare nested models (k = 1, 2, 3 components) and outputs comprehensive diagnostic plots, parameter estimates, and summary tables for downstream interpretation.

---

## Input data format

The pipeline expects a **tab‑separated file** containing one row per introgressed segment. A minimal example (`input/demo_Den_tracts.CHM13.txt`) is provided in this repository.

### Required columns (must be present with exact names):

| Column name           | Description                                                  |
| --------------------- | ------------------------------------------------------------ |
| `Assembly_ID`         | Unique identifier of the haploid assembly (e.g., `C021-CHA-S01_2-Pat`). |
| `weighted_match_rate` | Match rate to the Altai Denisovan reference (continuous between 0 and 1). |
| `width_seg`           | Length of the introgressed segment (in base pairs). Used for filtering. |
| `posterior`           | Posterior probability of introgression (e.g., from ASMAI). Used for filtering. |

Any other columns (e.g., `#chro`, `start`, `end`) are ignored by the GMM function but can be retained for reference.

---

## The GMM function: `fun_getGMM_result`

The core analysis is implemented in a single R function `fun_getGMM_result(df, Pop_para)`.

### Parameters:

- `df` : a **data frame** containing the required columns (see above) for **one haploid assembly** (i.e., all rows belong to the same `Assembly_ID`).
- `Pop_para` : a character string used as a label for output file naming (usually the `Assembly_ID`).

The function **does not return a value** but writes multiple output files to a specified directory (default: `~/Documents/00/`). It also generates diagnostic plots saved to the same location.

### Workflow overview

1. **Data filtering (implicitly)** – before calling the function, the user should apply desired filters (e.g., `weighted_match_rate > 0.3`, `width_seg > 20000`, `posterior >= 0.9`). The example notebook shows typical filters.
2. **k = 1 model**: computed manually (mean, sd) and BIC stored.
3. **k = 2 and k = 3 models**: fitted with `mixtools::normalmixEM` (parameters sorted by mean). BIC and log‑likelihood are extracted.
4. **Likelihood‑ratio tests (LRT)**:  
   - k=2 vs k=1: LRT statistic = 2 × (logLik₂ – logLik₁) ~ χ²(df = 3).  
   - k=3 vs k=2: similarly.
5. **Diagnostic plots**:
   - Histogram + fitted density curves for k = 1, 2, 3 (individual PDFs and a combined plot).
   - BIC values across k (line plot with ΔBIC and LRT p‑values annotated).
   - Bar plot showing mixing weights (λ) per component with μ and σ labels.
6. **Output files** (all saved to `~/Documents/00/`):
   - `{Pop_para}_denisovan_GMManalysis_Fit{k}.pdf` – fitted curves for each k.
   - `{Pop_para}_combined_gmm_fits_denisovan.pdf` – all fits stacked.
   - `{Pop_para}_GMM_BIC_mixtools.pdf` – BIC line plot.
   - `{Pop_para}_GMM_BIC_LRT_sig_test.txt` – table with BIC, ΔBIC, LRT p‑values per k.
   - `{Pop_para}_GMM_k1BIC_LRT_sig_test.txt` – BIC for k=1 (used later for merging).
   - `{Pop_para}_GMM_meanSD.txt` – estimated parameters (μ, σ, λ) for each k.
   - `{Pop_para}_GMM_lambda_plot.png` – bar plot of component proportions.
   - `{Pop_para}_MR_cluster.txt` – original data with assigned cluster labels (optional, commented out).

---

## How to run the analysis

### Step 1: Prepare the input data

Load your introgression tract file (e.g., from ASMAI) and apply quality filters. For example:

```r
library(tidyverse)

infile <- "/path/to/your/intro_tracts.bed.gz"
df_asmai <- read_tsv(infile) %>%
  mutate(width_seg = end - start, posterior = intro_prob)

# Join with population metadata if desired
df_popInfo <- read_tsv("population_info.txt")
df0 <- df_asmai %>% left_join(df_popInfo, by = "Assembly_ID")

df <- df0 %>%
  mutate(Denisovan = weighted_match_rate) %>%
  filter(Denisovan > 0.30,
         width_seg > 20000,
         posterior >= 0.9)
```

### Step 2: Call the function for each assembly

The function is designed to be applied per assembly. Use `group_walk()` or a loop:

```r
source("01.Gaussian_mixture models_function.Rmd")  # or source the script containing the function

df %>%
  group_by(Assembly_ID) %>%
  group_walk(~ {
    tryCatch({
      fun_getGMM_result(.x, .y$Assembly_ID)
    }, error = function(e) {
      message("Error in ", .y$Assembly_ID, ": ", e$message)
    })
  })
```

All output files will be written to `~/Documents/00/`. Make sure this directory exists or modify the output path inside the function.

---

## Output organisation and post‑processing

After running the function for many assemblies, the output folder (`~/Documents/00/`) will be cluttered with many files. A helper script (provided in the notebook) moves them into structured subfolders:

```r
# Create folders
folders <- c("plot1_bic", "plot2_gmm", "plot3_lambda", "plot4_bic_mclust",
             "plot5_density", "res1_bic", "res2_meanSD", "res3_optimal_bic_mclust",
             "res4bic_k1")
walk(folders, ~ dir_create(path("~/Documents/00", .x)))

# Move files based on suffix
# (see the notebook for the full mapping)
```

After organisation, two summary tables are created by merging individual results:

- `summary_bic/001_matchRate_BIC.txt` – contains BIC and LRT p‑values for all assemblies and all k.
- `summary_meanSD/001_GMM_meanSD.txt` – contains μ, σ, λ for all assemblies and all k.

These merged tables are used to determine the **optimal number of components** per assembly.

---

## Determining the optimal k (number of components)

The optimal model is selected using a **hierarchical likelihood‑ratio test**:

1. Compare k=3 vs k=2: if p < 0.05 → k=3 is accepted.
2. Else, compare k=2 vs k=1: if p < 0.05 → k=2 is accepted.
3. Otherwise, k=1 is chosen.

This logic is implemented in the notebook (section `step5: Calculate optimal BIC – directly based on p-value judgment`). It takes the merged BIC table and adds a column `k_fit` with the chosen number of components.

```r
df_bicMerge <- read_tsv("summary_bic/001_matchRate_BIC.txt")
df_optimal <- df_bicMerge %>%
  group_by(Pop) %>%   # Pop here corresponds to Assembly_ID
  summarise(...)       # (see the notebook for full details)
```

The result is a table with the best‑fit model for each assembly, ready for geographic or population‑level visualisation.

---

## Notes and recommendations

- **Minimum sample size**: The function includes checks; it will skip assemblies with too few segments (e.g., if after filtering fewer than ~50 segments remain). Adjust the threshold inside the function if needed.
- **Convergence issues**: `mixtools::normalmixEM` may occasionally fail to converge. The function handles these cases gracefully by returning `NA` for the problematic k and continuing.
- **Customising output directory**: The output path is hardcoded as `~/Documents/00/`. Change this by modifying the `str_c` calls inside `fun_getGMM_result`.
- **Parallelisation**: For many assemblies, consider using `furrr` or `parallel` to speed up the loop.

---

## Citation

If you use this pipeline, please cite our manuscript (**Deciphering complete archaic introgression sequences in modern human genomes**) and the underlying R packages:

- Benaglia, T., Chauveau, D., Hunter, D. R., & Young, D. S. (2009). mixtools: An R package for analyzing finite mixture models. *Journal of Statistical Software*, 32(6), 1–29.
