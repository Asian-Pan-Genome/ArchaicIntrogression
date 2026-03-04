import msprime
from math import exp
from matplotlib import pyplot as plt
import demesdraw
import os

'''
	Step1 : Simulate Demography
'''

seed = 1234     # Random Seed
out_dir = "Sim_Result"
if not os.path.exists(out_dir):
	os.makedirs(out_dir)
sample_dir = out_dir + "/sample"
if not os.path.exists(sample_dir):
	os.makedirs(sample_dir)

n_nonAFR = 100   # num. sampled non-African haplotypes
n_Nean = 2      # num. sampled Neand. haplotypes
n_AFR = 40      # num. sampled African haplotypes
n_bp = 100_000_000 # simulated 100 Mb
recomb_rate = 1.0e-8    # recombination rate
mu_rate = 2.4e-8        # mutation rate of base simulation 


N_NEAN = 1500 # Neanderthal effective size
N_AFR = 15000 # Ancestral African eff. size
N_OOA = 2000  # Out-of-African eff. seize
growth_rate = 0.02 # growth rate since agricultureA

T_NEAN_HUMAN = 16000 # neaderthal-human split time
T_OOA = 2400 # time of Out of Africa event
T_INTRO_START = 2000 # start of introgression
T_INTRO_END = T_INTRO_START - 20 # end of introgression
T_growth = 200 # time of recent growth

N_nonAFR_now = N_OOA * exp(T_growth * growth_rate)
N_AFR_now = N_AFR * exp(T_growth * growth_rate)

# migration rates are proportion of population made of new immigrants each generation
# Note: only one generation!!
# Real proportion = migration rate * migration time
M_AFR_nonAFR = 1e-10 # migration rate between AFR and non-AFR
M_INTRO = 0.0015 # Introgression Rate from Neanderthal into non-AFR 

print("\n>>> Set Demography ...")

demography = msprime.Demography()
demography.add_population(
	name="nonAFR",
	description="non-African(EUR or EAS) current pop.size",
	initial_size=N_nonAFR_now,
	growth_rate=growth_rate
)
demography.add_population(
	name="Neanderthal",
	description="Neanderthal pop.size",
	initial_size=N_NEAN
)
demography.add_population(
	name="AFR",
	description="African current pop.size",
	initial_size=N_nonAFR_now,
	growth_rate=growth_rate
)


# # Current Migration Rate Between AFR and non-AFR
demography.set_symmetric_migration_rate(["AFR", "nonAFR"], M_AFR_nonAFR)

# onset of growth with agriculture
demography.add_population_parameters_change(
	time=T_growth,
	initial_size=N_AFR,
	population="AFR",
	growth_rate=0
)
demography.add_population_parameters_change(
	time=T_growth,
	initial_size=N_OOA,
	population="nonAFR",
	growth_rate=0
)

# End of Introgression 【Note the direction!!】
demography.add_migration_rate_change(
	time=T_INTRO_END,
	rate=M_INTRO,
	source="nonAFR",
	dest="Neanderthal"
)

# Start of Introgression
demography.add_migration_rate_change(
	time=T_INTRO_START,
	rate=0,
	source="nonAFR",
	dest="Neanderthal"
)

# Human out of Africa
demography.add_mass_migration(
	time=T_OOA,
	source="nonAFR",
	dest="AFR",
	proportion=1.0
)
demography.add_symmetric_migration_rate_change(
	time=T_OOA,
	populations=["AFR", "nonAFR"],
	rate=0
) # Note: the migration rate needs to clear!!!

# Neanderthal split from Africa
demography.add_mass_migration(
	time=T_NEAN_HUMAN,
	source="Neanderthal",
	dest="AFR",
	proportion=1.0
)

#demography.debug()

# Plot the demography

graph = msprime.Demography.to_demes(demography)
fig, ax = plt.subplots()
demesdraw.tubes(graph, ax=ax, log_time=True, seed=seed)
plt.savefig(f"{out_dir}/note.deme.png")

fig, ax = plt.subplots()
demesdraw.size_history(graph, invert_x=True, ax=ax, log_size=True)
plt.savefig(f"{out_dir}/note.deme_size.png")


'''
	Step2 : Simulate Ancestry & Add Mutation
'''

print("\n>>> Simulate Ancestry ...")

# Define the sample list
samples = [
	msprime.SampleSet(n_nonAFR,population="nonAFR",ploidy=1),
	msprime.SampleSet(n_Nean,population="Neanderthal",ploidy=1,time=T_INTRO_START),
	msprime.SampleSet(n_AFR,population="AFR",ploidy=1),
]

# Run simulation
ts = msprime.sim_ancestry(
	samples=samples,
	demography=demography,
	sequence_length=n_bp,
	recombination_rate=recomb_rate,
	random_seed=seed
)

#print(ts.tables.nodes)

# Add mutation
ts = msprime.sim_mutations(
	ts,
	rate=mu_rate,           # mutation rate
	model=msprime.JC69()    # mutation model
)

print("\n>>> Get Simulate VCF file ...")

# Output vcf file
total_vcf = f"{out_dir}/note.sim.total.vcf"
with open(total_vcf,"w") as vcf_file:
	ts.write_vcf(vcf_file)

'''
	Step3 : Get Introgressed Region in bed format (truth dataset)
'''

def node_get_pop(tree, node, admix_time, split_time):
	while tree.get_time(node) <= admix_time:
		node = tree.get_parent(node)
	if tree.get_time(node) > split_time:
		return -1
	else:
		return tree.get_population(node)

def merge_adjacent_intervals(sorted_intervals):
	if not sorted_intervals:
		return []
	
	merged = [list(sorted_intervals[0])]
	
	for current in sorted_intervals[1:]:
		last = merged[-1]
		if current[0] == last[1]:
			# Adjacent region, merge
			last[1] = current[1]
		else:
			merged.append(list(current))
			
	return [tuple(interval) for interval in merged]

# key: nonAFR Sample ID
# value: introgressed interval
sample_intervals = {sample_id: [] for sample_id in range(n_nonAFR)}
nonAFRhaps = range(n_nonAFR)

print("\n>>> Get Introgressed Region ...")

# transverse all the tree
for tree in ts.trees():
	start = int(tree.interval.left)
	end = int(tree.interval.right)
	
	# Detect for the current region, which sample was introgressed
	intro_samples = []
	for sample_id in range(n_nonAFR):
		if node_get_pop(tree, sample_id, T_INTRO_START, T_NEAN_HUMAN) == 1:
			sample_intervals[sample_id].append((start, end))

# Write the introgressed-bed for every sample
for sample_id in range(n_nonAFR):
	# Merge Region
	merged = merge_adjacent_intervals(sample_intervals[sample_id])
	
	filename = f"{sample_dir}/nonAFR_{sample_id}.intro"
	
	with open(filename, "w") as f:
		f.write("#CHROM\tSTART\tEND\n")
		for interval in merged:
			f.write(f"1\t{interval[0]}\t{interval[1]}\n")
			
	print(f">>> nonAFR Sample {sample_id}: {len(merged)} merged intervals written")


'''
	Step4 : Get modified VCF for ASMaid HMM input
'''

print("\n>>> Get modified VCF for ASMaid HMM input ...")

with open(total_vcf,"r") as f:
	header_lines = []
	for line in f:
		if line.startswith("#CHROM"):
			headers = line.strip().split("\t")
			headers[9:] = ["Neanderthal"] + \
							 [f"African_{i+1}" for i in range(int(n_AFR/2))]
			break
		header_lines.append(line)
	
	# Store the remaining lines
	data_lines = list(f)

	for i in range(n_nonAFR):
		output_file = open(f"{sample_dir}/nonAFR_{i}.vcf","w")
		output_file.write("".join(header_lines))
		output_file.write("\t".join(headers) + "\n")

		for line in data_lines:
			if line.startswith("#"):
				continue
				
			fields = line.strip().split("\t")
			chrom, pos, id_, ref, alt, qual, filt, info, format_ = fields[:9]
			
			# Only contain bi-allele SNP
			if "," in alt:
				continue
				
			nonAFR_gt = fields[9+i]
			nean_afr_gt = fields[9+n_nonAFR:9+n_nonAFR+n_Nean+n_AFR]
			
			if nonAFR_gt == "1":
				# Flip the allele in Nean and AFR
				nean_afr_gt = [str(1-int(x)) for x in nean_afr_gt]
				# Flip the ref and alt
				tmp = ref
				ref = alt
				alt = tmp
			
			if all(x == "0" for x in nean_afr_gt):
				# Skip the line if all Nean and AFR are 0
				continue
				
			update_gt = []
			for j in range(int((n_Nean+n_AFR)/2)):
				gt = nean_afr_gt[2*j] + "/" + nean_afr_gt[2*j+1]
				update_gt.append(gt)
			
			new_fields = [chrom, pos, id_, ref, alt, qual, filt, info, format_] + update_gt
			output_file.write("\t".join(new_fields)+"\n")

		output_file.close()
		print(f">>> nonAFR Sample {i}: VCF modified !")
