#!/usr/bin/sh


sim_dir="Sim_Result"
archaic_vcf=$sim_dir/note.sim.total.diploid.nean.vcf
modern_vcf=$sim_dir/note.sim.total.diploid.modern_human.vcf
Result_dir=$sim_dir/Result_IBDmix

date
# Step 1 Generate Genotype
generate_gt -a $archaic_vcf -m $modern_vcf -o $Result_dir/genotype.vcf

# Step 2 Run IBDmix
ibdmix -g $Result_dir/genotype.vcf -s $sim_dir/sample.nonAFR.list -n Neanderthal -o $Result_dir/ibd_output.txt 
# Step 3 Filter Result
cat $Result_dir/ibd_output.txt | awk '$4-$3+1>50000 && $5>4' > $Result_dir/ibd_output.filtered.txt

echo "Finish!"
date

