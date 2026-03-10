import sys
import pandas as pd 
import numpy as np 

sample_order_file = sys.argv[1]
all_detail_file = sys.argv[2]
out_file = sys.argv[3]

sample_order_list = []
with open(sample_order_file, "r") as f:
	for line in f:
		sample = line.strip().split("\t")[0]
		sample_order_list.append(sample)

df = pd.read_csv(all_detail_file, sep="\t", header=None, comment='#', usecols = [3,4], names = ['block_id', 'sample'])

# construct sample x block matrix
matrix = pd.crosstab(df['sample'], df['block_id']).clip(upper=1)
# reindex to ensure all samples are included in the specified order
matrix = matrix.reindex(sample_order_list).fillna(0)

# print(matrix)

results = []

# Record for current samples, each block appears how many times
# index : block_id
current_counts = pd.Series(0, index=matrix.columns)

for i, sample_name in enumerate(sample_order_list, 1):
	print(f"Processing sample {i}/{len(sample_order_list)}: {sample_name}")
	current_counts += matrix.loc[sample_name]

	# Get active blocks (appear at least once in current samples)
	active_blocks = current_counts[current_counts > 0]
	total_active_blocks = len(active_blocks)

	frequencies = active_blocks / i 

	def get_label(count, freq):
		if count == 1 : return "singleton"
		if freq < 0.01 : return "rare"
		if freq < 0.05 : return "low_freq"
		if freq < 0.4 : return "common"
		return "high_freq"

	labels = pd.Series([get_label(c,f) for c,f in zip(active_blocks,frequencies)]).value_counts()

	res_entry = {
		'SampleNum' : i,
		'ChunkNum' : total_active_blocks,
		'Singleton' : labels.get('singleton', 0),
		'Rare' : labels.get('rare', 0),
		'LowFreq' : labels.get('low_freq', 0),
		'Common' : labels.get('common', 0),
		'HighFreq' : labels.get('high_freq', 0),
	}

	results.append(res_entry)

result_df = pd.DataFrame(results)
result_df.to_csv(out_file, sep="\t", index=False)

print("Done!")