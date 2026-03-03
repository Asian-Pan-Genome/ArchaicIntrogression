#!/bin/bash
chrom=`cat demo.info | cut -f1`
st=`cat demo.info | cut -f2`
mk=`cat demo.info | cut -f3`
echo ${chrom}
#awk 'OFS="\t" {print $1,$2-100000,$2+100000}' demo.info > demo.bed
#python3 ~/pyscript/sel_cal.py -ehh -chr ${chrom} -db ../../${chrom}/sv.${chrom}.db -snp -SV -indel -up 242696752 -l list -sam sample.list > ehh.sv.info
#cat <(awk '$9 > 0.05' ehh.sv.info | cut -f1-5,10 | sort | uniq | awk 'OFS="\t" {print $1,$2,$3,$4,$5,$2,$3,$5,1,$6}') ehh.sv.info | sort -k3,3 -k6,6n | awk '$9 > 0.05' > ehh.draw.info

region=`awk '{print $1":"$2-1000000"-"$2+1000000}' demo.info`
bcftools view -r ${region} -S <(grep EAS /share/home/zhanglab/user/nielei/project/APG/haplotype/hap_selection/sv_ehh/case6/sample.info | cut -f1) --force-samples /share/home/project/zhanglab/Primate_gene_loss/hg38_variants/1kg_3202/00.vcf/02.Biallelic_SNP_phased_panel/1kGP_high_coverage_Illumina.${chrom}.filtered.SNV_phased_panel.final.vcf.gz | bcftools norm -m+any | bcftools view -m 2 -M 2 -O z -o ehh.vcf.gz
Rscript ../ehh.R --marker ${mk}
mv Rplots.pdf ${chrom}_${st}.EAS.pdf
region=`awk '{print $1":"$2-1000000"-"$2+1000000}' demo.info`
bcftools view -r ${region} -S <(grep CHB /share/home/zhanglab/user/nielei/project/APG/haplotype/hap_selection/sv_ehh/case6/sample.info | cut -f1) --force-samples /share/home/project/zhanglab/Primate_gene_loss/hg38_variants/1kg_3202/00.vcf/02.Biallelic_SNP_phased_panel/1kGP_high_coverage_Illumina.${chrom}.filtered.SNV_phased_panel.final.vcf.gz | bcftools norm -m+any | bcftools view -m 2 -M 2 -O z -o ehh.vcf.gz
Rscript ../ehh.R --marker ${mk}
mv Rplots.pdf ${chrom}_${st}.CHB.pdf
region=`awk '{print $1":"$2-1000000"-"$2+1000000}' demo.info`
bcftools view -r ${region} -S <(grep CHS /share/home/zhanglab/user/nielei/project/APG/haplotype/hap_selection/sv_ehh/case6/sample.info | cut -f1) --force-samples /share/home/project/zhanglab/Primate_gene_loss/hg38_variants/1kg_3202/00.vcf/02.Biallelic_SNP_phased_panel/1kGP_high_coverage_Illumina.${chrom}.filtered.SNV_phased_panel.final.vcf.gz | bcftools norm -m+any | bcftools view -m 2 -M 2 -O z -o ehh.vcf.gz
Rscript ../ehh.R --marker ${mk}
mv Rplots.pdf ${chrom}_${st}.CHS.pdf
#cat <(bcftools view -v snps ehh.vcf.gz) <(bcftools view -H -r chr1:168829064 ehh.vcf.gz) | bcftools sort -O z -o ehh.snp.vcf.gz

#source /share/home/zhanglab/user/liujing/miniconda3/bin/activate /share/home/zhanglab/user/liujing/miniconda3/envs/repeatmasker
#RepeatMasker -pa 4  -species "Homo sapiens" -nolow -e ncbi -dir ./ /share/home/zhanglab/user/liuanguo/1.haplo_block/07.pos_sel/11.xp-nSL/slc19a2.insertion.pm100bp.fa
