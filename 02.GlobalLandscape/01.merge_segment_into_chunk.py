import sys
from collections import defaultdict

block_id = 0

def merge_intervals(intervals):
	global block_id
	''' merge overlapping or contiguous intervals '''
	if not intervals:
		return []
	sorted_intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
	merged = []
	start, end = sorted_intervals[0]
	for i in range(1, len(sorted_intervals)):
		s, e = sorted_intervals[i]
		if s <= end: # Overlapping or contiguous intervals
			end = max(end, e)
		else:
			block_id += 1
			merged.append((start, end, block_id))
			start, end = s, e
	block_id += 1
	merged.append((start, end, block_id))  # Add last interval
	return merged


if len(sys.argv) != 4:
	print(f"Usage: python {sys.argv[0]} <bed_list_file> <probability_fiter> <out_prefix>")
	sys.exit(1)
bed_list_file = sys.argv[1]
# bed_list_file: A file containing sample names and their corresponding BED files, formatted as:	
# <sample_name>\t<bed_file_path>
prob_filt = float(sys.argv[2])
out_prefix = sys.argv[3] 
# Output prefix for the merged results
# --- two files will be generated:
# 1. {out_prefix}.merged.bed: Merged intervals with weighted match rates & sample numbers
# 2. {out_prefix}.merged.detail: Detailed statistics for each merged region

sample_bed_list = []

sample_data = defaultdict(list)  # {sample_name: [(chr, start, end, ... , match-rate), ...]}
all_intervals = defaultdict(list)  # store all segments for merging
sample_num = 0

with open(bed_list_file,"r") as f:
	for line in f:
		sample_num += 1
		sample_name, bed_file = line.strip().split('\t')
		sample_bed_list.append((sample_name,bed_file))

# Read each sample's BED file and store data
for sample_name, bed_file in sample_bed_list:
	with open(bed_file, 'r') as f:
		for line in f:
			if line.startswith("#"):
				continue
			fields = line.strip().split('\t')
			prob = float(fields[4])
			# Filter by probability
			if prob <= prob_filt: 
				continue
			# Store Tracts
			chrom, start, end, seg_id, gt_intro_num, cnv_intro_num = fields[0], int(fields[1]), int(fields[2]), fields[3], fields[5], fields[6]
			match_rate = float(fields[7])
			sample_data[sample_name].append((chrom, start, end, seg_id, gt_intro_num, cnv_intro_num, match_rate))
			all_intervals[chrom].append((start, end))

# Generate merged intervals
merged_regions = []  # Format: [(chr, start, end)]
for chrom, intervals in all_intervals.items():
	for start, end, merge_id in merge_intervals(intervals):
		merged_regions.append((chrom, start, end, merge_id))

# Calculate overlap for each sample
overlap_records = []  # Store final output records
merged_regions_info = []  # Store merged region info for output

for m_chr, m_start, m_end, m_id in merged_regions:
	total_overlap = 0
	sample_details = []  # Store details for current merged region
	
	# Calculate total overlap length for the merged region
	for sample, intervals in sample_data.items():
		for (chrom, start, end, seg_id, gt_intro_num, cnv_intro_num, match_rate) in intervals:
			if chrom != m_chr: 
				continue
			overlap_len = max(0, min(end, m_end) - max(start, m_start))
			if overlap_len > 0:
				total_overlap += overlap_len
				sample_details.append({
					"sample": sample,
					"start": start,
					"end": end,
					"seg_id": seg_id,
					"gt_intro_num": gt_intro_num,
					"cnv_intro_num": cnv_intro_num,
					"match_rate": match_rate,
					"overlap": overlap_len
				})
	# Calculate weighted match rate for the merged region
	weighted_values = []
	sample_list_set = set()
	for detail in sample_details:
		sample_list_set.add(detail["sample"])
		weight = detail["overlap"] / total_overlap
		weighted_val = weight * detail["match_rate"]
		weighted_values.append(weighted_val)
		# Add record to output list (with metadata)
		overlap_records.append((
			m_chr, m_start, m_end, m_id,
			detail["sample"],
			detail["start"], detail["end"], detail["seg_id"],
			detail["gt_intro_num"],detail["cnv_intro_num"],
			detail["match_rate"],
			weight,weighted_val
		))
	sample_list_num = len(sample_list_set)
	sample_list_freq = sample_list_num / sample_num
	weighted_match_rate = sum(weighted_values)
	# Append merged region record
	merged_regions_info.append((m_chr, m_start, m_end, m_id, weighted_match_rate, sample_list_num, sample_list_freq))

# Write merged results to output files
merged_bed_file = f"{out_prefix}.merged.bed"
merged_detail_file = f"{out_prefix}.merged.detail"
with open(merged_bed_file, 'w') as bed_out, open(merged_detail_file, 'w') as detail_out:
	bed_out.write("#chr\tstart\tend\tblock_id\tweighted_match_rate\tsample_list_num\tsample_list_freq\n")
	for m_chr, m_start, m_end, m_id, weighted_match_rate, sample_list_num, sample_list_freq in merged_regions_info:
		bed_out.write(f"{m_chr}\t{m_start}\t{m_end}\t{m_id}\t{weighted_match_rate:.4f}\t{sample_list_num}\t{sample_list_freq:.4f}\n")
	
	detail_out.write("#chr\tstart\tend\tblock_id\tsample\tsample_start\tsample_end\tsample_seg_id\tsample_gt_intro_num\tsample_cnv_intro_num\tsample_match_rate\tsample_weight\tsample_weighted_match_rate\n")
	for m_chr, m_start, m_end, m_id, sample, s_start, s_end, s_id, s_gt, s_cnv, s_match_rate, s_weight, s_weight_match_rate in overlap_records:
		detail_out.write(f"{m_chr}\t{m_start}\t{m_end}\t{m_id}\t{sample}\t{s_start}\t{s_end}\t{s_id}\t{s_gt}\t{s_cnv}\t{s_match_rate:.4f}\t{s_weight:.4f}\t{s_weight_match_rate:.4f}\n")