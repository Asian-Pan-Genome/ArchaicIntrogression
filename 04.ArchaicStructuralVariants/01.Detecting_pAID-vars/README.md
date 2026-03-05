**Detailed descriptions are provided in the manuscript.**

# Requirements
Before running the pipeline, you should make sure these software/packages are installed:
- bcftools
- snakemake
- vcfwave
- pysam
- truvari
- pandas
- numpy
- scipy
- BWA
- samtools
- 


# Pangenome graph construction for each introgression chunk
@mingyu

Questions: In which step is the `detail` file generated?? This file would be used in the sec `Detecting pAID-vars` (see below).



After constructing the graph, there would be the `$chunk.vcf.gz` file containing variants decomposed from the graph, which could be used as input for the subsequent pipeline.

# Variants decomposition
Here, the `$chunk.vcf.gz` file decomposed from the graph above was subsequently processed using the [VCF preparation pipeline](https://github.com/eblerjana/genotyping-pipelines/tree/main/prepare-vcf-MC) and a [SV collapsing pipeline](https://github.com/Han-Cao/collapse-bubble).

First, copy the folder [scripts](https://github.com/Asian-Pan-Genome/ArchaicIntrogression/tree/main/04.ArchaicStructuralVariants/01.Detecting_pAID-vars/scripts) in your work directory, and edit the [config file](https://github.com/Asian-Pan-Genome/ArchaicIntrogression/blob/main/04.ArchaicStructuralVariants/01.Detecting_pAID-vars/scripts/config.yaml) in the `scripts` folder, as well as provide a TSV file specifying sex (1=male, 2=female) of each sample.

Next, you can make the `variants_decomposition` pipeline for each chunk:
```
bash scripts/variants_decompose.sh $home $threads $list $REF
# then, run bash/shell script in each chunk folder
```

Finally, you will find `${chunk}_filtered_ids_biallelic.sort.uniqid.vcfwave.sort.merge.sort.dedup.vcf.gz` file in each chunk folder, which obtains the de-redundancy variants set.


# Inversions calling
Although `vcfwave` in the above pipe could identify a few inversions, we additionally integrated inversion calls from [PAV](https://github.com/EichlerLab/pav) and [LSGvar](https://github.com/Hanjunmin/LSGvar).

@feifei

This process would generate the file `INV.tsv` containing all inversion calls.
# Detecting pAID-vars
By integrating inversion calls from the above step, you can now identify pAID-vars by comparing the allele frequencies between introgressed and non-introgressed haplotypes using Fisher's exact test (or Chi-square test), for each variant.
```
python scripts/pAID-vars_freq_test.py $input.detail $chunk INV.tsv > $chunk.freq_test.tsv
```
After completing, you would get a table like this:
```
chrom   pos     block_id        chi_square_pvalue       fisher_pvalue   introgressed_0  introgressed_1  introgressed_missing    non_introgressed_0      non_introgressed_1      non_introgressed_missing        SVTYPE  SVLEN   REF     ALT
chr5    79323784        1616    0.10191297541610678     0.0474567241451096      138     0       0       459     13      0                       G       A
chr5    79323853        1616    5.272859186855923e-06   1.1496582132236e-07     136     2       0       392     80      0       INS     2       C       CCA
chr5    79323853        1616    0.10191297541610678     0.0474567241451096      138     0       0       459     13      0       INS     4       C       CCACA
chr5    79323853        1616    1.0     1.0     138     0       0       470     2       0       INS     6       C       CCACACA
chr5    79324126        1616    0.014136843028468286    0.002447463715531566    138     0       0       448     24      0                       T       C
chr5    79324199        1616    1.0     1.0     138     0       0       471     1       0                       T       C
chr5    79324226        1616    0.512544559202569       0.22622950819672136     137     1       0       472     0       0                       C       T
chr5    79324269        1616    0.8047633853003324      1.0     138     0       0       469     3       0                       C       T
```
where `chi_square_pvalue` and `fisher_pvalue` are unadjusted p-values, and you can use some thresholds to filter significant ones, that is, pAID-vars. Here, in the manuscript, we select the threshold $5e-8$, which is the common choice for GWAS research.

Then, you could concat all tables into archaic-level `$archaic.freq_test.tsv`.


# Genotyping pAID-vars in Neanderthal and Denisovan genomes
Here we defined the pAID-vars present in at least one Neanderthal or Denisovan genomes were highly-confident (HC). For small variants (pAID-SMVs, including SNVs and InDels), we could map archaic reads to the reference genome and use GATK to call variants, while for SVs, which are too difficult to accuartely calling, validation was performed by analyzing mapping coverage profiles and sequence clipping signals.

## Archaic reads mapping and variants calling
`BWA-ALN`
@mingyu

This step would generate file `$archaic.vcf.gz` for downstream analysis.

## Validation for pAID-SNVs
```
python scripts/validation_pAID-SNVs.py $archaic.freq_test.tsv $archaic.vcf.gz > $archaic.snv.confirmed.tsv
```

## Validation for pAID-InDels
For the same indel variants, different callers and methods would generate inconsistent representations (similar to SVs). Here we useed `truvari bench` for comparision.
### Classifying pAID-InDels based on the allele frequency
Since the reference genome may represent either a modern or an archaic allele at specific loci, we should employ different validation strategies. Therefore, here we classify them based on the allele frequency present in both introgressed and non-introgressed haplotypes, with the results from the script `pAID-vars_freq_test.py`.
```
python scripts/classify_pAID-InDels.py $archaic.freq_test.tsv $archaic.freq_test.tsv.A.small_ins.vcf $archaic.freq_test.tsv.A.small_del.vcf $archaic.freq_test.tsv.B.small_ins.vcf $archaic.freq_test.tsv.B.small_del.vcf
```

### Comparing variant alleles using `truvari bench`
For class A, where the reference represents a non-archaic allele, we first construct a `pseudo-vcf ($archaic.fix_gt.vcf.gz)` by fixing genotypes for archaic genome to `1/1` for all variant records. You could easily create one using `awk` or `sed`:
```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Nean-Altai
chr1    12567   chr1_12567_C_A  C       A       415     .       AF=0.5;AQ=415;AN=6;AC=2 GT:DP:AD:GQ:PL:RNC      1/1:10:0,10:28:415,30,0:..
chr1    23463   chr1_23463_T_C  T       C       109     .       AF=0.333333;AQ=109;AN=4;AC=2    GT:DP:AD:GQ:PL:RNC      1/1:10:0,10:28:415,30,0:..
chr1    52084   chr1_52084_A_G  A       G       18      .       AF=0.125;AQ=18;AN=6;AC=1        GT:DP:AD:GQ:PL:RNC      1/1:10:0,10:28:415,30,0:..
```
For class B, where the reference represents an archaic allele, we construct `$archaic.hom_alt.vcf.gz`, only extracting variants with homozygous alternative genotypes (`1/1`) called from at least one archaic genome. You could easily create one using `bcftools`:
```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Nean-Altai
chr1    12567   chr1_12567_C_A  C       A       415     .       AF=0.5;AQ=415;AN=6;AC=2 GT:DP:AD:GQ:PL:RNC      1/1:10:0,10:28:415,30,0:..
chr1    23463   chr1_23463_T_C  T       C       109     .       AF=0.333333;AQ=109;AN=4;AC=2    GT:DP:AD:GQ:PL:RNC      1/1:10:0,10:28:415,30,0:..
chr1    52362   chr1_52362_A_T  A       T       295     .       AF=0.333333;AQ=295;AN=4;AC=2    GT:DP:AD:GQ:PL:RNC      1/1:10:0,10:28:415,30,0:..
```

Then, truvari was used for comparing the identified pAID-InDels and variants called based on archaic reads.
```
for i in ins del; do
    rm -rf $archaic.freq_test.tsv.A.small_${i}/
    truvari bench -b $archaic.fix_gt.vcf.gz -c $archaic.freq_test.tsv.A.small_${i}.vcf.gz -o $archaic.freq_test.tsv.A.small_${i}/ -f $REF -r 250 -p 0.5 -P 0.5 -t --pick multi -d -s 1 -S 1
    truvari refine -t 16 -f $REF $archaic.freq_test.tsv.A.small_${i}/

    rm -rf $archaic.freq_test.tsv.B.small_${i}/
    truvari bench -b $archaic.hom_alt.vcf.gz -c $archaic.freq_test.tsv.B.small_${i}.vcf.gz -o $archaic.freq_test.tsv.B.small_${i}/ -f $REF -r 250 -p 0.5 -P 0.5 -t --pick multi -d -s 1 -S 1
    truvari refine -t 16 -f $REF $archaic.freq_test.tsv.B.small_${i}/
done
```

### Validation
```
python scripts/validation_pAID-InDels.py $archaic.freq_test.tsv $archaic.freq_test.tsv.A.small_ins/refine.comp.vcf.gz $archaic.freq_test.tsv.A.small_del/refine.comp.vcf.gz $archaic.freq_test.tsv.B.small_ins/refine.comp.vcf.gz $archaic.freq_test.tsv.B.small_del/refine.comp.vcf.gz $archaic.indel.confirmed.tsv
```

## Validation for pAID-SVs
### Archaic reads mapping
While `BWA-ALN` could not generate clipping signals which are important for SVs identification, we additionally utilized `BWA-MEM` for archiaic reads mapping. 
@mingyu

### Validation
First, collect mapping depth information from `$archaic.BWA-ALN.bam`. Before validation, please run `samtools depth` to get the whole-genome avarge depth `${archaic.BWA-ALN.depth}` for archaic genome.

Then, validate pAID-SVs using only depth profile. 
```
python scripts/validation_pAID-SVs.BWA-ALN.py $archaic.freq_test.tsv $archaic.BWA-ALN.bam 500000 32 ${archaic.BWA-ALN.depth} $archaic.BWA-ALN # please specify your window size, here is 500Kb
```




For structural variants (pAID-SVs), validation was performed by analyzing mapping coverage profiles and sequence clipping signals from archaic genomes. Archaic WGS data were aligned to T2T-CHM13 using BWA-ALN and BWA-MEM (v0.7.17) for collecting mapping depth and clipping information, respectively. When T2T-CHM13 represents a non-archaic allele, a pAID-SV was defined as HC if the mapping depth was less than 50% of the whole-genome average or if clipping signals were observed. Due to the low coverage of the Vindija and Chagyrskaya genomes, only the Altai Neanderthal was utilized for depth-based assessments. When T2T-CHM13 represents an archaic allele, a pAID-SV was classified as HC if the mapping depth exceeded 50% of the whole-genome average or if the count of clipping reads were less than half of the whole-genome depth.


