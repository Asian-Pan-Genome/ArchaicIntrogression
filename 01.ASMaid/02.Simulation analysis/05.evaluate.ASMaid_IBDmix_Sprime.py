def cal_evaluation_matrix(true_bed_list,true_len,call_file):
	overlap_len = 0
	call_len = 0

	with open(call_file,"r") as f:
		for line in f:
			if line.startswith("#"):
				continue
			fields = line.strip().split("\t")
			start, end = int(fields[1]), int(fields[2])
			call_len += end - start + 1
			for t_start, t_end in true_bed_list:
				if min(end,t_end) - max(start,t_start) > 0:
					# Overlap
					overlap_len += min(end,t_end) - max(start,t_start)
	# print(overlap_len)
	precision_len = overlap_len / call_len
	recall_len = overlap_len / true_len
	f1_len = (2 * precision_len * recall_len) / (precision_len + recall_len)
	false_positive_rate = (true_len - overlap_len) / 100000000
	return precision_len, recall_len, f1_len, false_positive_rate


Sim_dir = "Sim_Result"
out = "05.eval.result"
sample_num = 100

with open(out,"w") as out_f:
	out_f.write(f"sample\tprecision\trecall\tf1\tfpr\tsoftware\n")
	for i in range(sample_num):
		print(f">>> Process non-afr {i}...")
		intro_file = f"{Sim_dir}/sample/nonAFR_{i}.intro"
		ASMaid_file = f"{Sim_dir}/Result_ASMaid/sample/nonAFR_{i}.ASMaid.afr10.prob09.bed"
		IBDmix_file = f"{Sim_dir}/Result_IBDmix/sample/nonAFR_{i}.ibdmix.filt.bed"
		Sprime_file = f"{Sim_dir}/Result_Sprime/sample/nonAFR_{i}.sprime.bed"
		true_bed_list = []
		true_len = 0
		with open(intro_file,"r") as true_f:
			for line in true_f:
				if line.startswith("#"):
					continue
				fields = line.strip().split("\t")
				start, end = int(fields[1]), int(fields[2])
				true_len += end - start + 1
				true_bed_list.append((start,end))
		# print(true_bed_list)
		pl, rl, fl, fpr = cal_evaluation_matrix(true_bed_list,true_len,ASMaid_file)
		out_f.write(f"nonAFR_{i}\t{pl:.4f}\t{rl:.4f}\t{fl:.4f}\t{fpr:.6f}\tASMaid\n")
		# print(pl, rl, fl, fpr)
		pl, rl, fl, fpr = cal_evaluation_matrix(true_bed_list,true_len,IBDmix_file)
		out_f.write(f"nonAFR_{i}\t{pl:.4f}\t{rl:.4f}\t{fl:.4f}\t{fpr:.6f}\tIBDmix\n")

		pl, rl, fl, fpr = cal_evaluation_matrix(true_bed_list,true_len,Sprime_file)
		out_f.write(f"nonAFR_{i}\t{pl:.4f}\t{rl:.4f}\t{fl:.4f}\t{fpr:.6f}\tSprime\n")

print("Finish All!")