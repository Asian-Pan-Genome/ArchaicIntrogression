home=$1
threads=$2
list=$3
REF=$4

cat $list | while read -r chunk; do
    cd $home/${chunk}
    
    echo "
        set -x
        set -oe pipefail
        python $home/scripts/fix_gap_alleles.py ${chunk}.vcf.gz ${chunk}.fix.vcf.gz
        mv ${chunk}.fix.vcf.gz ${chunk}.vcf.gz
        bcftools index -f -t --threads $threads ${chunk}.vcf.gz
        
        rm -rf .snakemake/ && rm -rf vcf/
        sed \"s/Den_1/${chunk}/g\" $home/scripts/config.yaml > config.yaml
        ln -sf $home/scripts/prepare-vcf-MC/workflow/
        snakemake -s $home/scripts/prepare-vcf-MC/workflow/Snakefile -j $threads
        bcftools sort --write-index -Oz -o ${chunk}_filtered_ids.sort.vcf.gz vcf/${chunk}/${chunk}_filtered_ids.vcf
        bcftools norm --threads $threads -f $REF -Ov vcf/${chunk}/${chunk}_filtered_ids_biallelic.vcf | bcftools sort --write-index -Oz -o ${chunk}_filtered_ids_biallelic.sort.vcf.gz
        
        python $home/scripts/annotate_var_id.py -i ${chunk}_filtered_ids_biallelic.sort.vcf.gz -o ${chunk}_filtered_ids_biallelic.sort.uniqid.vcf.gz
        bcftools annotate -x INFO/AT -Ov ${chunk}_filtered_ids_biallelic.sort.uniqid.vcf.gz | bcftools +fill-tags -- -t AC,AN,AF | vcfwave -t $threads -I 1000 | bgzip -@ $threads -c > ${chunk}_filtered_ids_biallelic.sort.uniqid.vcfwave.vcf.gz
        bcftools norm --threads $threads -f $REF ${chunk}_filtered_ids_biallelic.sort.uniqid.vcfwave.vcf.gz | bcftools sort --write-index -Oz -o ${chunk}_filtered_ids_biallelic.sort.uniqid.vcfwave.sort.vcf.gz

        python $home/scripts/collapse_bubble.py -i ${chunk}_filtered_ids_biallelic.sort.uniqid.vcfwave.sort.vcf.gz -o ${chunk}_filtered_ids_biallelic.sort.uniqid.vcfwave.sort.merge.vcf.gz --map ${chunk}_filtered_ids_biallelic.sort.uniqid.vcfwave.sort.merge.mapping.txt
        bcftools sort --write-index -Oz -o ${chunk}_filtered_ids_biallelic.sort.uniqid.vcfwave.sort.merge.sort.vcf.gz ${chunk}_filtered_ids_biallelic.sort.uniqid.vcfwave.sort.merge.vcf.gz

        python $home/scripts/merge_duplicates.py -i ${chunk}_filtered_ids_biallelic.sort.uniqid.vcfwave.sort.merge.sort.vcf.gz -o ${chunk}_filtered_ids_biallelic.sort.uniqid.vcfwave.sort.merge.sort.dedup.vcf.gz -c repeat
        bcftools +fill-tags ${chunk}_filtered_ids_biallelic.sort.uniqid.vcfwave.sort.merge.sort.dedup.vcf.gz -Oz -o ${chunk}_filtered_ids_biallelic.sort.uniqid.vcfwave.sort.merge.sort.dedup.1.vcf.gz -- -t AC,AN,AF 
        mv ${chunk}_filtered_ids_biallelic.sort.uniqid.vcfwave.sort.merge.sort.dedup.1.vcf.gz ${chunk}_filtered_ids_biallelic.sort.uniqid.vcfwave.sort.merge.sort.dedup.vcf.gz
        python $home/scripts/fix_info.py ${chunk}_filtered_ids_biallelic.sort.uniqid.vcfwave.sort.merge.sort.dedup.vcf.gz ${chunk}_filtered_ids_biallelic.sort.uniqid.vcfwave.sort.merge.sort.dedup.fix.vcf.gz
        mv ${chunk}_filtered_ids_biallelic.sort.uniqid.vcfwave.sort.merge.sort.dedup.fix.vcf.gz ${chunk}_filtered_ids_biallelic.sort.uniqid.vcfwave.sort.merge.sort.dedup.vcf.gz
        
        rm -rf .snakemake/ && rm -rf vcf/
    " > ${chunk}.sh
done
