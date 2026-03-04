
# cM = bp × recombination_rate
sim_dir = "Sim_Result/"
vcf_file = sim_dir + "note.sim.total.diploid.vcf"
genetic_map_file = sim_dir + "note.sim.genetic_map.txt"

recombination_rate = 1.0e-8 
positions = []
with open(vcf_file) as f, open(genetic_map_file, "w") as out_f:
	for line in f:
		if not line.startswith("#"):
			chrom, pos, *_ = line.split()
			cM = int(pos) * recombination_rate * 100
			out_f.write(f"1 rs{pos} {cM:.6f} {pos}\n")
