# ASMaid detects complete archaic introgression sequences from haplotype-resovled assemblies
This folder includes workflow and code scripts for the article section 1 (ASMaid detects complete archaic introgression sequences from haplotype-resovled assemblies). As for ASMaid pipeline, this folder contains the mapping and variant calling steps to work as pre-processing for the main pipeline. **The main ASMaid pipeline is available at** https://github.com/Asian-Pan-Genome/ASMaid. Simulation analysis and CN-based case analysis are also included in this folder.

## 01. Mapping and variant calling
This pipeline handles DNA alignment and variant discovery for both modern HGDP (Bergström et al., 2020, Science) samples and archaic hominins (Neanderthal Altai, Neanderthal Chagyrskaya, and Denisovan). Ten east African samples from HGDP were randomly selected for the introgression background control to mitigate ILS.  
To ensure high-quality variant calls across different types of sequencing data:
 - **Archaic DNA**: Processed via `blw aln` (optimized for damaged reads) and called via `GATK HaplotypeCaller`.
 - **Modern DNA**: Processed via `bwa mem` and called via `DeepVariant`.
 - **Joint Genotyping**: All samples are merged using `Glnexus` (DeepVariantWGS config).

## 02. Simulation analysis

This module is dedicated to generating simulation data, validating the performance of our proposed method (**ASMaid**) and comparing it against state-of-the-art introgression detection tools (**IBDmix** and **Sprime**). 

### Key Components:
- **Demographic Simulation**: Utilizing `msprime` to simulate complex Out-of-Africa evolutionary scenarios, including Neanderthal-to-Human gene flow (archaic introgression).
- **Cross-Tool Benchmarking**: A unified pipeline to process simulated data through multiple detection algorithms.
- **Statistical Evaluation**: Robust calculation of Precision, Recall, F1-score, and False Positive Rate (FPR) by comparing inferred tracts with the simulation ground truth.

## 03. CN-based case analysis (_MGAM_ gene)

## Key Tools & Versions