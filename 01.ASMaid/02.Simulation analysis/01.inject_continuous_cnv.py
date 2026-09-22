import bisect
import sys
import numpy as np
import pandas as pd

# ----------------------------------------------------------------------
# CLI Arguments and Path Parsing
# ----------------------------------------------------------------------
if len(sys.argv) < 4 or len(sys.argv) > 5:
	print(
		"Usage: python inject_continuous_cnv_custom.py <input.onlyGT.vcf>"
		" <input_gt.bed> <out_prefix> [--plot]"
	)
	sys.exit(1)

CHROM_MAX_LEN = 100_000_000

vcf_in = sys.argv[1]
gt_bed_in = sys.argv[2]
out_prefix = sys.argv[3]

vcf_out = f"{out_prefix}.addDP.vcf"
bed_out = f"{out_prefix}.addDP.intro"
detailed_bed_out = f"{out_prefix}.addDP.intro.detail"

do_plot = False
if len(sys.argv) == 5:
	flag = sys.argv[4].lower()
	if flag in ["--plot", "-plot", "plot", "true"]:
		do_plot = True
		plot_out = f"{out_prefix}.tract_plot.pdf"
	else:
		print(f"[Warning] Unrecognized plotting flag '{sys.argv[4]}', skipping PDF generation.")

ARCHAIC_IDX = 0
AFR_START = 1
AFR_NUM = 20

np.random.seed(42)


# ----------------------------------------------------------------------
# 0. Depth (DP) Generation Function Based on Copy-Number Zygosity
# ----------------------------------------------------------------------
def generate_cnv_dp(cnv_type, is_archaic_sample=False):
	"""Generates continuous normalized Read Depth (DP) modeled on copy-number states and Gaussian noise:

	- Normal (CN=2):           DP ~ Normal(1.0, std)
	- Heterozygous DEL (CN=1): DP ~ Normal(0.48, std) [50% probability]
	- Homozygous DEL (CN=0):   DP ~ Normal(0.02, 0.02) [50% probability, background sequencing noise]
	- Duplication (CN=3/4/5):  DP ~ Normal(1.6, std) [60%], Normal(2.2, std) [30%], or Normal(2.8, std) [10%]
	"""
	std = 0.25 if is_archaic_sample else 0.12

	if cnv_type == "del":
		# 50% Homozygous Deletion (CN=0, mean=0.02), 50% Heterozygous Deletion (CN=1, mean=0.48)
		is_hom = np.random.rand() < 0.6
		if is_hom:
			dp = np.random.normal(0.02, 0.02)
		else:
			dp = np.random.normal(0.48, std)

	elif cnv_type == "dup":
		# 60% Single Copy Duplication (CN=3, mean=1.6), 40% Multi-Copy Duplication (CN=4/5, mean=2.2/2.8)
		rand_val = np.random.rand()
		if rand_val < 0.6:
			mean_val = 1.6
		elif rand_val < 0.9:
			mean_val = 2.2
		else:
			mean_val = 2.8
		dp = np.random.normal(mean_val, std)

	else:
		# Normal state (CN=2, mean=1.0)
		dp = np.random.normal(1.00, std)

	return max(0.0, dp)


# ----------------------------------------------------------------------
# 1. Construct Archaic CNV Regions and Background Population Noise
# ----------------------------------------------------------------------
gt_df = pd.read_csv(
	gt_bed_in,
	sep="\t",
	header=None,
	comment="#",
	dtype={0: str},
)
gt_df = gt_df.iloc[:, :3]
gt_df.columns = ["chrom", "start", "end"]

cnv_regions = []
detailed_records = []

# (1) Record original archaic introgressed Ground Truth tracts
for _, row in gt_df.iterrows():
	detailed_records.append(
		(str(row["chrom"]), int(row["start"]), int(row["end"]), "GT_Original")
	)

# (2) Inject archaic-derived CNV signals into ~40% of GT introgressed tracts
for _, row in gt_df.iterrows():
	chrom, start, end = str(row["chrom"]), int(row["start"]), int(row["end"])
	seg_len = end - start

	if np.random.rand() < 0.4:
		# 1. CNV length sampling: bounded within [1,000 bp, 40,000 bp] proportional to seg_len
		ratio = np.random.uniform(0.3, 1.5)
		cnv_len = int(np.clip(seg_len * ratio, 1000, 40000))
		
		# 2. Centroid offset: CNV midpoint deviates from GT midpoint via Gaussian distribution (std = 0.3 * seg_len)
		gt_center = (start + end) / 2.0
		shift = np.random.normal(0, seg_len * 0.3)
		cnv_center = gt_center + shift
		
		# 3. Derive physical start and end positions from the continuous probability model
		cnv_start = max(0, int(cnv_center - cnv_len / 2.0))
		cnv_end = max(cnv_start + 1000, int(cnv_center + cnv_len / 2.0))

		cnv_type = np.random.choice(["del", "dup"])
		cnv_regions.append({
			"chrom": chrom,
			"start": cnv_start,
			"end": cnv_end,
			"type": cnv_type,
			"is_archaic": True,
			"carrier_mask": None,
		})
		detailed_records.append(
			(chrom, cnv_start, cnv_end, f"CNV_Archaic_{cnv_type}")
		)

# (3) Append 50 background polymorphic CNV noise regions (ILS / Non-introgressed Population CNV)
chrom_name = str(gt_df.iloc[0]["chrom"])

for _ in range(50):
	rand_start = int(np.random.uniform(100, CHROM_MAX_LEN - 30000))
	rand_end = rand_start + int(np.random.uniform(500, 5000))
	cnv_type = np.random.choice(["del", "dup"])

	# Pre-assign population carrier masks (Allele Frequency: 15% - 40%)
	carrier_freq = np.random.uniform(0.15, 0.40)
	carrier_mask = np.random.rand(500) < carrier_freq

	cnv_regions.append({
		"chrom": chrom_name,
		"start": rand_start,
		"end": rand_end,
		"type": cnv_type,
		"is_archaic": False,
		"carrier_mask": carrier_mask,
	})
	detailed_records.append(
		(chrom_name, rand_start, rand_end, f"CNV_Noise_{cnv_type}")
	)

# ---- Merge overlapping intervals to build the Combined Ground Truth BED ----
all_intervals = []
for _, row in gt_df.iterrows():
	all_intervals.append(
		(str(row["chrom"]), int(row["start"]), int(row["end"]))
	)
for region in cnv_regions:
	if region["is_archaic"]:
		all_intervals.append((region["chrom"], region["start"], region["end"]))

all_intervals.sort(key=lambda x: (x[0], x[1]))

merged_gt = []
for interval in all_intervals:
	if not merged_gt:
		merged_gt.append(interval)
	else:
		prev_chrom, prev_start, prev_end = merged_gt[-1]
		curr_chrom, curr_start, curr_end = interval
		if curr_chrom == prev_chrom and curr_start <= prev_end:
			merged_gt[-1] = (prev_chrom, prev_start, max(prev_end, curr_end))
		else:
			merged_gt.append(interval)

with open(bed_out, "w") as f_bed:
	f_bed.write("#chrom\tstart\tend\n")
	for chrom, m_start, m_end in merged_gt:
		f_bed.write(f"{chrom}\t{m_start}\t{m_end}\n")
		detailed_records.append((chrom, m_start, m_end, "GroundTruth_Merged"))

with open(detailed_bed_out, "w") as f_det:
	f_det.write("#chrom\tstart\tend\ttype\n")
	for chrom, s, e, t in sorted(detailed_records, key=lambda x: (x[0], x[1])):
		f_det.write(f"{chrom}\t{s}\t{e}\t{t}\n")

print(f">>> Combined BED generated: {bed_out}")
print(f">>> Detailed Records BED generated: {detailed_bed_out}")

# ----------------------------------------------------------------------
# 2. Construct O(log N) Binary Search Index for Quick Interval Lookup
# ----------------------------------------------------------------------
cnv_regions.sort(key=lambda x: x["start"])
cnv_starts = [c["start"] for c in cnv_regions]


def get_cnv_info_fast(pos):
	idx = bisect.bisect_right(cnv_starts, pos) - 1
	if idx >= 0:
		region = cnv_regions[idx]
		if region["start"] <= pos <= region["end"]:
			return region
	return None


# ----------------------------------------------------------------------
# 3. Process Input VCF Line-by-Line & Inject Simulated DP Information
# ----------------------------------------------------------------------
print(f">>> Injecting read depth (DP) fields into VCF...")

with open(vcf_in, "r") as fin, open(vcf_out, "w") as fout:
	for line in fin:
		if line.startswith("##FORMAT=<ID=GT"):
			fout.write(line)
			fout.write(
				'##FORMAT=<ID=DP,Number=1,Type=Float,Description="Normalized'
				' Read Depth">\n'
			)
			continue
		if line.startswith("#"):
			fout.write(line)
			continue

		fields = line.strip().split("\t")
		pos = int(fields[1])
		fields[8] = "GT:DP"

		cnv_info = get_cnv_info_fast(pos)
		samples = fields[9:]
		new_samples = []

		for idx, s_val in enumerate(samples):
			gt_val = s_val.split(":")[0]
			is_archaic_sample = idx == ARCHAIC_IDX

			if cnv_info is None:
				# Normal genomic regions: no CNV events
				dp = generate_cnv_dp(
					"normal", is_archaic_sample=is_archaic_sample
				)

			else:
				cnv_type = cnv_info["type"]
				is_archaic_cnv = cnv_info["is_archaic"]

				if is_archaic_cnv:
					# ==================================================
					# Case A: Archaic Introgressed CNVs
					# ==================================================
					if is_archaic_sample:
						# Archaic reference genome baseline (Normal Depth ~ 1.0)
						dp = generate_cnv_dp(
							"normal", is_archaic_sample=is_archaic_sample
						)

					elif AFR_START <= idx < (AFR_START + AFR_NUM):
						# Non-introgressed outgroup (AFR): lacks archaic haplotype, displays DEL/DUP state
						dp = generate_cnv_dp(
							cnv_type, is_archaic_sample=is_archaic_sample
						)

					else:
						# Modern human target samples: 85% carry introgressed segment (Normal DP), 15% non-carriers (Variant DP)
						is_inherited = np.random.rand() < 0.85
						if is_inherited:
							dp = generate_cnv_dp(
								"normal", is_archaic_sample=is_archaic_sample
							)
						else:
							dp = generate_cnv_dp(
								cnv_type, is_archaic_sample=is_archaic_sample
							)

				else:
					# ==================================================
					# Case B: Background Polymorphic CNVs (ILS / Population Noise)
					# ==================================================
					is_carrier = cnv_info["carrier_mask"][idx]
					if is_carrier:
						dp = generate_cnv_dp(
							cnv_type, is_archaic_sample=is_archaic_sample
						)
					else:
						dp = generate_cnv_dp(
							"normal", is_archaic_sample=is_archaic_sample
						)

			new_samples.append(f"{gt_val}:{dp:.2f}")

		fout.write("\t".join(fields[:9] + new_samples) + "\n")

print(f">>> VCF depth injection completed: {vcf_out}")

# ----------------------------------------------------------------------
# 4. Optional: Generate Visualization PDF for Genomic Tracts
# ----------------------------------------------------------------------
if do_plot:
	print(f">>> Generating genomic tract visualization plot: {plot_out} ...")
	import matplotlib.pyplot as plt

	track_mapping = {
		"GT_Original": (1, "#377EB8", "Original GT"),
		"CNV_Archaic_del": (2, "#E41A1C", "Archaic CNV (DEL)"),
		"CNV_Archaic_dup": (2, "#FF7F00", "Archaic CNV (DUP)"),
		"CNV_Noise_del": (3, "#999999", "Noise CNV (DEL)"),
		"CNV_Noise_dup": (3, "#CCCCCC", "Noise CNV (DUP)"),
		"GroundTruth_Merged": (4, "#4DAF4A", "Combined GT"),
	}

	fig, ax = plt.subplots(figsize=(14, 4.5), dpi=300)
	seen_labels = set()

	for chrom, s, e, t in detailed_records:
		if t in track_mapping:
			y_idx, color, label_name = track_mapping[t]
			start_mb = s / 1e6
			end_mb = e / 1e6
			width = max(end_mb - start_mb, 0.02)

			lbl = label_name if label_name not in seen_labels else ""
			seen_labels.add(label_name)

			ax.barh(
				y=y_idx,
				width=width,
				left=start_mb,
				height=0.55,
				color=color,
				edgecolor="none",
				alpha=0.85,
				label=lbl,
			)

	ax.set_yticks([1, 2, 3, 4])
	ax.set_yticklabels(
		["Original GT", "Archaic CNV", "Noise CNV", "Combined GT"],
		fontsize=10,
		fontweight="bold",
	)
	ax.set_xlabel("Genomic Position (Mbp)", fontsize=11, fontweight="bold")
	ax.set_title(
		"Genomic Tract Distribution across 0-100 Mbp",
		fontsize=13,
		fontweight="bold",
	)
	ax.set_xlim(0, CHROM_MAX_LEN / 1e6)
	ax.grid(axis="x", linestyle="--", alpha=0.5)

	handles, labels = ax.get_legend_handles_labels()
	ax.legend(
		handles,
		labels,
		loc="upper right",
		bbox_to_anchor=(1.0, 1.25),
		ncol=3,
		frameon=True,
	)

	plt.tight_layout()
	plt.savefig(plot_out, format="pdf", bbox_inches="tight")
	plt.close()
	print(f">>> Genomic tract plot successfully saved to: {plot_out}")