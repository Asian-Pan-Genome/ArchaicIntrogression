#!/usr/bin/sh

sim_dir="Sim_Result"
# Step 0 Prepare Genetic Map
python 04.call_Sprime_for_simulation.00_get_genetic_map.py

# Step 1 Run Sprime
sprime="java -jar /share/home/zhanglab/user/suomingyu/src/software/sprime.jar"
$sprime gt=$sim_dir/note.sim.total.diploid.modern_human.vcf \
    outgroup=$sim_dir/sample.AFR.list \
    map=$sim_dir/note.sim.genetic_map.txt \
    out=$sim_dir/Result_Sprime/sprime_output_mu24 \
    mu=2.4e-8

# Step 2 Get Match & Mismatch Info
echo "Get Match & Mismatch Info"
python 04.call_Sprime_for_simulation.02_get_match_info.py \
    $sim_dir/Result_Sprime/sprime_output_mu24.score \
    $sim_dir/note.sim.total.diploid.nean.vcf

echo "Finish!"



