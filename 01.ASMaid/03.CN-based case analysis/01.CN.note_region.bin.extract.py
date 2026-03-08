import sys
import math

if __name__ == "__main__":
	# input
	bin_num = int(sys.argv[1]) # 100 -> 100 bins [by base pairs]
	region = sys.argv[2]  # chr:start-end
	vcf_pas = sys.argv[3] # eg: /share/home/zhanglab/user/suomingyu/Archaic/HGDP-Archaic-APG/Result/C001-CHA-E01-Mat/GATK
	sample = sys.argv[4] # assembly_name, use to store in the 1st column of the output file && the output file name
	Archaic_index = int(sys.argv[5]) # 0->Nean-Altai,1->Nean-Cha,2->Den-Den
	Archaic_lst = ["Nean-Altai","Nean-Cha","Den-Den"]

	# Store the mean norm_DP of each bin
	# if the bin is not covered by any CNV, the mean_norm_DP is 1
	chro = region.split(":")[0]
	region_start = int(region.split(":")[1].split("-")[0])
	region_end = int(region.split(":")[1].split("-")[1])
	bin_size = math.ceil((region_end - region_start) / bin_num)
	bin_norm_DP_Archaic = [[] for _ in range(bin_num)]
	bin_norm_DP_AFR = [[] for _ in range(bin_num)]


	vcf = f"{vcf_pas}/{chro}/Merge.{chro}.sample.mask.bi-allele.normDP.vcf"
	with open(vcf, "r") as f:
		for line in f:
			if line.startswith("#"):
				continue
			parts = line.strip('\n').split('\t')
			pos = int(parts[1])
			if pos < region_start or pos > region_end:
				continue
			bin_index = (pos - region_start) // bin_size
			bin_norm_DP_Archaic[bin_index].append(float(parts[Archaic_index+9].split(':')[-1]))
			# AFR sample : from column 12 -> 12+10, get avarge norm_DP
			bin_norm_DP_AFR[bin_index].append(sum([float(parts[i+12].split(':')[-1]) for i in range(10)]) / 10)
	#print(bin_norm_DP_Archaic)

	# Calculate the mean norm_DP of each bin
	for i in range(bin_num):
		bin_norm_DP_Archaic[i] = sum(bin_norm_DP_Archaic[i]) / len(bin_norm_DP_Archaic[i]) if bin_norm_DP_Archaic[i] else 1
		bin_norm_DP_AFR[i] = sum(bin_norm_DP_AFR[i]) / len(bin_norm_DP_AFR[i]) if bin_norm_DP_AFR[i] else 1
	# Write the results to the output file
	output_file = f"{sample}.{Archaic_lst[Archaic_index]}.{bin_num}bin.result.txt"
	with open(output_file, "w") as f:
		# sample type(Archaic or AFR)  bin1  bin2  bin3  ...
		f.write("#sample\ttype\t"+"\t".join([f"bin{i+1}" for i in range(bin_num)])+"\n")
		f.write(f"{sample}\tArchaic\t" + "\t".join(map(lambda x: f"{x:.2f}", bin_norm_DP_Archaic[:bin_num])) + "\n")
		f.write(f"{sample}\tAFR\t" + "\t".join(map(lambda x: f"{x:.2f}", bin_norm_DP_AFR[:bin_num])) + "\n")

