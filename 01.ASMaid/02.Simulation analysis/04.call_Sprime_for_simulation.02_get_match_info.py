import sys

sprime_output = sys.argv[1]
archaic_vcf = sys.argv[2]
sprime_output_add_match = sprime_output + ".add_match_info"

archaic_snp_dict = dict()

with open(archaic_vcf,"r") as f:
	for line in f:
		if line.startswith("#"):
			continue
		fields = line.strip().split("\t")
		pos = fields[1]
		gt = fields[9].split("|")
		archaic_snp_dict[pos] = gt

with open(sprime_output,"r") as infile, open(sprime_output_add_match,"w") as outfile:
	for line in infile:
		if line.startswith("CHROM"):
			line = line.strip() + "\tN_MATCH\n"
		else:
			fields = line.strip().split("\t")
			pos = fields[1]
			archaic_allele = fields[6]
			if archaic_allele in archaic_snp_dict[pos]:
				line = line.strip() + "\tmatch\n"
			else:
				line = line.strip() + "\tmismatch\n"
		outfile.write(line)

