import rich.progress
import os

out_dir = "Sim_Result"
total_vcf = f"{out_dir}/note.sim.total.vcf"
out_vcf = f"{out_dir}/note.sim.total.diploid.vcf"

nean_vcf = f"{out_dir}/note.sim.total.diploid.nean.vcf"
modern_human_vcf = f"{out_dir}/note.sim.total.diploid.modern_human.vcf"

# Method A: Hap1 + Hap2 -> Diploid
# Method B: Hap + Hap -> Diploid 

n_Nean = 2
n_AFR = 40
n_nonAFR = 100

# Method A for Neanderthal and African samples
# Method B for non-African samples

with rich.progress.open(total_vcf,"r") as infile, open(out_vcf, "w") as outfile:
	for line in infile:
		if line.startswith("#"):
			if line.startswith("#CHROM"):
				headers = line.strip().split("\t")
				headers[9:] = [f"nonAFR_{i}" for i in range(n_nonAFR)] + \
							["Neanderthal"] + \
							[f"African_{i+1}" for i in range(n_AFR//2)]
				outfile.write("\t".join(headers) + "\n")
			else:
				outfile.write(line)
		else:
			fields = line.strip().split("\t")
			# Only contain bi-allele SNP
			if "," in fields[4]:
				continue
			# Combine Hap1 and Hap2 for each individual
			new_fields = fields[:9]
			new_fields += [f"{fields[9+i]}|{fields[9+i]}" for i in range(n_nonAFR)]
			new_fields += [f"{fields[9+n_nonAFR+2*i]}|{fields[9+n_nonAFR+2*i+1]}" for i in range((n_Nean+n_AFR)//2)] 
			outfile.write("\t".join(new_fields) + "\n")

# Split the Combined VCF into Archaic & Modern Human VCFs
os.system(f"bcftools view -s Neanderthal -Ov -o {nean_vcf} {out_vcf}")
os.system(f"bcftools view -s ^Neanderthal -Ov -o {modern_human_vcf} {out_vcf}")

# Get List of Non-African and AFR Samples
os.system(f"bcftools query -l {out_vcf} > {out_dir}/sample.list")
os.system(f"cat {out_dir}/sample.list | grep nonAFR > {out_dir}/sample.nonAFR.list")
os.system(f"cat {out_dir}/sample.list | grep African > {out_dir}/sample.AFR.list")