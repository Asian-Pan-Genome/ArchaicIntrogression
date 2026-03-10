# Global landscape of complete archaic introgression sequences
## Overview
This pipeline processes the raw `introgression segments` inferred by **ASMaid** to construct high-confidence `introgressed chunks`. The process involves merging overlapping segments, removing potential Incomplete Lineage Sorting (ILS) noise, and classifying the remaining chunks based on their population frequency. Finally, it evaluates the saturation of the archaic pangenome as a function of sample size.

## Data Resources
This repository provides high-confidence archaic introgression datasets, processed from 611 assemblies (from HPRCy1, HGSVC3, APGp1, SEA3K projects and one melanesian) using the T2T-CHM13 reference coordinate system. All segments have been filtered for high posterior probability (> 0.9) and cleaned of ancestral Incomplete Lineage Sorting (ILS).

### File Links & Descriptions
| File Name | Description |
| ---------- | ---------- |
| [Neanderthal.CHM13.prob09.rmILS.merged.bed](data/Neanderthal.CHM13.prob09.rmILS.merged.bed) | Merged Neanderthal-introgressed chunks (BED format).|
| [Neanderthal.CHM13.prob09.rmILS.merged.detail](data/Neanderthal.CHM13.prob09.rmILS.merged.detail) | Compositional details for each Neanderthal-introgressed chunk.|
| [Denisovan.CHM13.prob09.rmILS.merged.bed](data/Denisovan.CHM13.prob09.rmILS.merged.bed) | Merged Denisovan-introgressed chunks (BED format).|
| [Denisovan.CHM13.prob09.rmILS.merged.detail](data/Denisovan.CHM13.prob09.rmILS.merged.detail) | Compositional details for each Denisovan-introgressed chunk.|

### File Format Specifications
**1. BED Files (`*.merged.bed`)**: These files define the genomic boundaries and population-level prevalence of the introgressed chunks.
+ `#chr`, `start`, `end`: Genomic coordinates based on the **T2T-CHM13** assembly.
+ `block_id`: Unique identifier for each merged introgression chunk.
+ `weighted_match_rate`: The match rate of the chunk, calculated as the weighted average of the match rates of its constituent segments, weighted by their respective lengths.
+ `sample_list_num`: Total number of samples carrying this specific chunk.
+ `sample_list_freq`: Prevalence of the chunk across the 611 assemblies.
  
**2. Detail Files (`*.merged.detail`)**: These files provide the granular composition of each chunk, mapping individual segments back to their supporting evidence.
+ `block_id`: Links the segment to the primary `block_id` defined in the BED file.
+ `sample`: Identifier of the individual sample carrying the segment.
+ `sample_start`, `sample_end`: Genomic coordinates of the segment in the T2T-CHM13 reference.
+ `sample_seg_id`: The identifier for the original segment source.
+ `sample_gt_intro_num`: Count of genotype-based (GT) supportive sites within the segment.
+ `sample_cnv_intro_num`: Count of copy number-based (CN) supportive sites within the segment
+ `sample_match_rate`: The raw match rate of this individual segment.
+ `sample_weight`: The ratio of this segment's length to the total length of the introgressed chunk.
+ `sample_weighted_match_rate`: The individual segment's match rate multiplied by its weight.

## Code Structure
+ `01.merge_segment_into_chunk.py`: Merges multiple samples' raw segments into discrete, non-overlapping chunks and calculates chunk-level weighted match rates.
+ `02.remove_ILS_segment.py`: Filters out chunks likely resulting from ancestral ILS based on African population frequency and length criteria.
+ `03.cumulative_stat_for_chunk_number.py`: Performs population-level categorization (singleton, rare, etc.) and saturation growth curve analysis.

## Detailed Usage
### Step 1: Merging Segments into Chunks
Segments from different individuals are merged based on chromosomal overlap.
```shell
python 01.merge_segment_into_chunk.py <bed_list_file> <probability_filter> <out_prefix>
```
+ **Parameters:**
  + `bed_list_file`: A tab-seperated file containing the sample names and their corresponding bed files for merging (formatted as: `<sample_name>\t<bed_file_path>`). The bed file was genenrated from `ASMaid` pipeline and then lifted to `T2T-CHM13` coordinate, detailed format could be found at https://github.com/Asian-Pan-Genome/ASMaid.
  + `probability_filter`: Float value to threshold the posterior probability of the input segments (e.g., 0.9).
+ **Outputs:**
  +  `<out_prefix>.merged.bed`: Contains merged intervals, weighted match rates, and sample count/frequency statistics.
  +  `<out_prefix>.merged.detail`: Detailed breakdown of each merged chunk’s composition.
  
### Step 2: Removing ILS Segments
To differentiate archaic introgression from ILS, we apply rigorous filtering based on African sample distribution. A chunk is removed due to ILS if:
+ **African Ancestry Predominance:** The proportion of African introgressed individuals among all carriers of the chunk is ≥ 60%.
+ **Localized LWK Enrichment:** ≥ 3 individuals in the African LWK population (Luhya in Webuye, Kenya) carry the chunk, and the cumulative length of the introgressed segment within these individuals exceeds 50% of the total chunk length.

```bash
python 02.remove_ILS_segment.py <merged_bed> <merged_detail> <AFR_sample_list> <LWK_sample_list> <output_prefix>
# Parameters:
## <merged_bed> and <merged_detail> are generated by step 1.
## <AFR_sample_list> and <LWK_sample_list> are lists of samples from the African and LWK populations.
```

### Step 3: Frequency Classification & Saturation Analysis
For cumulative chunk counts, we utilized the final 610-sample aggregate as a reference set, determining the proportion of the total `“introgressed pangenome”` captured with each incremental addition.

+ **Frequency Categories:**
  + `Singleton`: introgressed sample number == 1
  + `Rare`: introgressed sample frequency (f) < 1%
  + `Low-frequency`: 1% ≤ f < 5%
  + `Common`: 5% ≤ f < 40%
  + `High-frequency`: f ≥ 40%

```bash
python 03.cumulative_stat_for_chunk_number.py <sample_order_file> <merged_detail> <output_file>
# Parameters:
## <sample_order_file> is a sample list file (seperated by <\n>) stores the order of the samples for cumulative analysis.
## <merged_detail> is generated by step 1.
```