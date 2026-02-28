# 01.Mapping and variant calling pipeline
## Overview
This section describes the automated pipeline for processing raw sequencing data (FASTQ) into a high-quality joint-genotyped cohort VCF. The pipeline accommodates two distinct data types: modern DNA (HGDP) and Archaic DNA (Neanderthal/Denisovan).

## Detailed Usage
### Step 1: Configure Data Root Directory & Input Files
The pipeline relies on a metadata list to batch-process samples. You have two options to provide this information:

#### Option A: Automatic Generation (Default)
The script `01.main_pipeline.sh` contains a function `generate_sample_config` that scans the `DATA_ROOT` to build a sample list.

- **HGDP List Format**: Files within `DATA_ROOT/HGDP.list/` (e.g., `Mbuti.list`) must follow this tab-delimited format:  
    `Population` [tab] `Sample_Name` [tab] `/Absolute/Path/To/Fastq_Folder`
- **Archaic Data**: FASTQ files for Neanderthal and Denisovan should be placed in `DATA_ROOT/Data/XXX`.
- **Note**: All FASTQ files are assumed to be **pre-processed/quality-controlled** (QC-passed).

#### Option B: Manual Input (Recommended for Custom Sets)
If you prefer to provide your own sample list or wish to bypass the automatic scanning:

1. **Comment out** the `generate_sample_config` line in `01.main_pipeline.sh`.
2. **Create a file** at `${CONFIG_DIR}/${ASSEMBLY}.sample_list.txt` with the following 4-column format:
    ```text
        #Format: Population  Name  Fastq_Path  Mode
        Mbuti  HGDP00450  /path/to/fastq_dir  mem
        Nean   Altai      /path/to/altai_dir  aln
    ```
   + **mem**: For paired-end modern DNA (using `bwa mem` + `DeepVariant`).
   + **aln**: For single-end archaic DNA (using `bwa aln` + `GATK`).


### Step 2: Mapping and Variant Calling
Execute the main controller. It uses a chromosome list generated at `${CONFIG_DIR}/${ASSEMBLY}.chro_list.txt` (extracted from the reference `.fai`) to parallelize calling.
```bash
# Usage: bash 00_main_pipeline.sh <Assembly_Name> <Ref_Fasta> <Project_Root> <Data_Root>
bash 00_main_pipeline.sh C115-CMG02-Mat /share/home/project/zhanglab/APG/Freezev0.9/C115-CMG02/Mat/C115-CMG02_Mat.v0.9.fasta ./project_dir /path/to/DATA_ROOT
```

### Step 3: Cohort Merging (Joint Genotyping)
After all individual jobs finish (signaled by `FinishDV` or `FinishGATK` flags), merge the results:
```bash
# Usage: bash 02_joint_genotype.sh <Assembly_Name> <Project_Root>
bash 02_joint_genotype.sh C115-CMG02-Mat ./project_dir
```
This step uses `GLnexus` to combine individual GVCFs into a single, analysis-ready cohort VCF for each chromosome.