# Pangenome graph construction for each introgression chunk
@mingyu

After constructing the graph, there would be the `$chunk.vcf.gz` file containing variants decomposed from the graph, which could be used as input for the subsequent pipeline.

# Variants decomposition
Here, the `$chunk.vcf.gz` file decomposed from the graph above was subsequently processed using the [VCF preparation pipeline](https://github.com/eblerjana/genotyping-pipelines/tree/main/prepare-vcf-MC) and a [SV collapsing pipeline](https://github.com/Han-Cao/collapse-bubble).

Before running the pipeline, you should make sure these softwares/packages installed:
- bcftools
- snakemake
- vcfwave
- pysam
- truvari
- pandas
- numpy

Then, copy the folder [scripts](https://github.com/Asian-Pan-Genome/ArchaicIntrogression/tree/main/04.ArchaicStructuralVariants/01.Detecting_pAID-vars/scripts) in your work directory, and edit the [config file](https://github.com/Asian-Pan-Genome/ArchaicIntrogression/blob/main/04.ArchaicStructuralVariants/01.Detecting_pAID-vars/scripts/config.yaml) in `scripts` folder, as well as provide a TSV file specifying sex (1=male, 2=female) of each sample.

After that, you can make the `variants_decomposition` pipeline for each chunk:
```
bash scripts/variants_decompose.sh $home $threads $list $REF
# then, run bash/shell script in each chunk folder
```

Finally, you will find `${chunk}_filtered_ids_biallelic.sort.uniqid.vcfwave.sort.merge.sort.dedup.vcf.gz` file in each chunk folder, which obtains the de-redundancy variants set.


# Inversions calling
Although `vcfwave` in the above pipe could identify a few inversions, we additionally integrated inversion calls from [PAV](https://github.com/EichlerLab/pav) and [LSGvar](https://github.com/Hanjunmin/LSGvar).

@feifei


# Detecting pAID-vars
