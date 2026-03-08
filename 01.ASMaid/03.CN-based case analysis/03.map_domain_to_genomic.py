import sys
import pandas as pd

if len(sys.argv) != 4:
	print(f"Usage: python {sys.argv[0]} <domtab_file> <gff_file> <output_file>")
	sys.exit(1)


# gff_file: use 02.result/*.tmp.gff

domtab_file = sys.argv[1]
gff_file = sys.argv[2]
output_file = sys.argv[3]

# 1. Get CDS coordinates from GFF file
cds_regions = []
strand = "+"
gene_start = 0

with open(gff_file, "r") as f:
	for line in f:
		if line.startswith("#"): continue
		parts = line.strip().split("\t")
		if parts[2] == "transcript" or parts[2] == "mRNA":
			gene_start, gene_end = int(parts[3]), int(parts[4])
			strand = parts[6]
		if parts[2] == "CDS":
			cds_regions.append((int(parts[3]), int(parts[4])))

if strand == "+":
	cds_regions.sort()
else:
	cds_regions.sort(reverse=True)

# 2. Decode domtab and convert the coordinate
results = []
with open(domtab_file, "r") as f:
	for line in f:
		if line.startswith("#") : continue
		p = line.split()
		domain_name, gene_id = p[0], p[3]
		aa_start, aa_end = int(p[17]), int(p[18]) # ali from/to

		# get genomic position (0-based)
		nt_start_offset = (aa_start - 1) * 3
		nt_end_offset = (aa_end * 3) - 1

		def get_genomic_pos(offset):
			current_len = 0
			for c_start, c_end in cds_regions:
				block_len = c_end - c_start + 1
				if current_len <= offset < current_len + block_len:
					if strand == "+":
						return c_start + (offset - current_len)
					else:
						return c_end - (offset - current_len)
				current_len += block_len
			return None
		g_start = get_genomic_pos(nt_start_offset)
		g_end = get_genomic_pos(nt_end_offset)

		if g_start and g_end:
			# let g_start be the smaller one, to plot easier
			real_start, real_end = min(g_start, g_end), max(g_start, g_end)
			rel_start = real_start - gene_start
			rel_end = real_end - gene_start
			results.append([gene_id, gene_start, gene_end, domain_name, aa_start, aa_end, real_start, real_end, rel_start, rel_end])

df = pd.DataFrame(results, columns=["gene_id", "gene_start", "gene_end", "domain_name", "aa_start", "aa_end", "g_start", "g_end", "rel_start", "rel_end"])
df.to_csv(output_file, sep="\t", index=False)
