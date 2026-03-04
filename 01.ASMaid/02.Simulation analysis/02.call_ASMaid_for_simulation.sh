#!/usr/bin/sh

sim_dir="/Sim_Result/sample"
ASMaid="python /share/home/zhanglab/user/suomingyu/Archaic/HGDP-Archaic-APG/Code/02.Admixture/assemblyAdmix-2/01.hmm.py"

for vcf in $(ls $sim_dir/*.vcf)
do
    outprefix="${vcf:0:-4}"
    echo "--------------------------------------"
    echo $outprefix
    # Test Different AFR Number [5,10,15,20]
    $ASMaid -i $vcf -o $outprefix.AFR-5 -t 32 --only-gt \
        -archaic_index 0 -afr_start 1 -afr_number 5
        
    $ASMaid -i $vcf -o $outprefix.AFR-10 -t 32 --only-gt \
        -archaic_index 0 -afr_start 1 -afr_number 10

    $ASMaid -i $vcf -o $outprefix.AFR-15 -t 32 --only-gt \
        -archaic_index 0 -afr_start 1 -afr_number 15
    
    $ASMaid -i $vcf -o $outprefix.AFR-20 -t 32 --only-gt \
        -archaic_index 0 -afr_start 1 -afr_number 20

    # Test Different Genotype Threshold [0,0.05,0.1,0.15]
    $ASMaid -i $vcf -o $outprefix.AFR-10.gt005 -t 32 --only-gt \
        -archaic_index 0 -afr_start 1 -afr_number 10 --gt-threshold 0.05
    
    $ASMaid -i $vcf -o $outprefix.AFR-10.gt010 -t 32 --only-gt \
        -archaic_index 0 -afr_start 1 -afr_number 10 --gt-threshold 0.1
    
    $ASMaid -i $vcf -o $outprefix.AFR-10.gt015 -t 32 --only-gt \
        -archaic_index 0 -afr_start 1 -afr_number 10 --gt-threshold 0.15    
done




