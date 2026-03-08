# 03. CN-based case analysis (_MGAM_ gene)
## Overview

This directory contains a suite of scripts for the analysis of the _MGAM_ gene locus. This locus is characterized by significant archaic introgression with multiple Copy Number (CN)-based supportive sites. Our analysis integrates bin-based depth profiling, protein domain mapping, and assembly-level synteny to resolve the structural landscape of _MGAM_ across different lineages.

## Code Structure
+ `01.CN.note_region.bin.extract.py`: Extracts normalized sequencing depth across specified bins for CNV visualization across different samples. In the article, 100 bins were set for this regions.
+ `02.if_Gap.predict_domain.sh`: Assesses assembly quality (Gap detection) and predicts Pfam domains using hmmscan in the target gene region (e. g. _MGAM_).
+ `03.map_domain_to_genomic.py`: Translates protein domain coordinates (AA) back to absolute genomic coordinates.
+ **Synteny Workflow**: (Commands included below) Uses `minimap2` and `SVbyEye` for structural alignment visualization.

## Detailed Usage
### Step 1: Bin-based Depth Profiling
To characterize the copy number states, we divide the introgression tract into equal-sized bins and extract the normalized depth (normDP).
```bash
# Usage: python 01.CN.note_region.bin.extract.py <bin_num> <region> <vcf_dir> <sample_id> <archaic_index>
## bin_num: 100 for practice in this article
## archaic_index: 0->Neanderthal Altai, 1->Neanderthal Chagyrskaya, 2->Denisovan
python 01.CN.note_region.bin.extract.py 100 chr7:1111000-1121000 ./VCF_Dir SampleA 0
```
#### Input VCF Format
The analysis requires a specific VCF format in `vcf_dir`, which also serves as the standard input for the **ASMaid** pipeline. The VCF must be structured such that archaic samples appear first, followed by modern human samples. Structure details are as follows:
  + **Sample Ordering**: Archaic samples must occupy the initial columns (e.g., `Nean-Altai`, `Nean-Cha`, `Den-Den`), followed by the modern human panel (e.g., HGDP samples).
  + **Metadata**: Each variant record includes a custom `NORM_DP` field within the `FORMAT` column, representing the normalized depth for that specific sample.
  + **Filtering**: The input VCFs are pre-filtered to include only high-quality bi-allelic SNPs within the target genomic regions.
    ```
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Nean-Altai	Nean-Cha	Den-Den	Mbuti-HGDP00983	Mbuti-HGDP00468	Mbuti-HGDP00450	BantuKenya-HGDP01408	BantuKenya-HGDP01406	BantuKenya-HGDP01413	San-HGDP00992	San-HGDP01029	BantuSouthAfrica-HGDP00993	BantuSouthAfrica-HGDP01031	Sardinian-HGDP01068	Sardinian-HGDP01072	Sardinian-HGDP00666	PapuanSepik-HGDP00546	PapuanSepik-HGDP00547	PapuanSepik-HGDP00542
    C001-CHA-E01#Mat#chr1	7261	C001-CHA-E01#Mat#chr1_7261_A_C	A	C	10	.	AF=0.052632;AQ=10;AC=18;AN=38	GT:DP:AD:GQ:PL:RNC:NORM_DP	0/0:0:0,0:0:0,0,0:..:0.00	0/1:0:0,0:0:0,0,0:..:0.00	0/0:0:0,0:0:0,0,0:..:0.00	0/1:0:0,0:0:0,0,0:..:0.00	0/1:2:2,0:6:0,6,59:..:0.08	1/1:2:0,2:2:10,14,0:..:0.08	0/1:0:0,0:0:0,0,0:..:0.00	0/1:0:0,0:0:0,0,0:..:0.00	0/1:0:0,0:0:0,0,0:..:0.00	0/1:0:0,0:0:0,0,0:..:0.00	0/1:0:0,0:0:0,0,0:..:0.00	0/1:0:0,0:0:0,0,0:..:0.00	0/1:1:1,0:3:0,3,29:..:0.04	0/1:0:0,0:0:0,0,0:..:0.00	0/1:0:0,0:0:0,0,0:..:0.00	0/1:0:0,0:0:0,0,0:..:0.00	0/1:27:27,0:81:0,81,809:..:1.08	0/1:21:21,0:33:0,33,569:..:0.81	0/1:24:24,0:72:0,72,719:..:0.96
    ```
Note: For details on how to generate these filtered, normalized VCFs from raw alignment data (BAMs), please refer to the pre-processing pipeline provided in the https://github.com/Asian-Pan-Genome/ASMaid.

### Step 2: Assembly Integrity & Functional Domain Annotation
We assess the continuity of the *MGAM* gene in the assemblies and predict its functional domains. Comparative analyses were restricted to assemblies containing the MANE Select transcript. Coding sequences (CDS) were extracted and translated via `gffread`, followed by protein domain identification using `HMMER` against the Pfam database ($E\text{-value} < 1 \times 10^{-5}$).
```bash
bash 02.if_Gap.predict_domain.sh input.gff assembly.fa output_prefix
# input.gff: gff file were annotated using Liftoff v1.6.3 from the GRCh38.p14 (GENCODE v47) reference
```

### Step 3: Genomic Coordinate Mapping
This step maps the identified HMMER protein domain hits back to the absolute genomic coordinates for precise visualization.
```shell
python 03.map_domain_to_genomic.py output.domtab temporary.gff mapping_result.txt
# output.domtab, temporary.gff were generated from 02.if_Gap.predict_domain.sh step
```

### Step 4: Synteny Analysis and Visualization
Compare different assemblies to visualize structural rearrangements like expansions or deletions (`Fig.1d` in article).

#### Alignment Generation
```bash
# Generate PAF file using minimap2
minimap2 -N 50 -p 0.1 -c -eqx $target.fa $query.fa > alignment.paf
```
#### Visualization (R - SVbyEye):
```R
# Plotting the alignment with SVbyEye
plotAVA(paf.table = paf,
        binsize = 1000,
        perc.identity.breaks = c(80, 90, 95, 99, 99.5, 99.9),
        seqnames.order = order)
```
