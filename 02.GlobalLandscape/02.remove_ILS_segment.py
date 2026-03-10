import sys
import pickle
import rich.progress

class MergedTract(object):

	def __init__(self):
		self.seg_len = 0
		self.sample_num = 0
		self.afr_num = 0
		self.lwk_num = 0
		self.lwk_len = 0
		self.temp_line_list = []

	def add_line(self,line,LWK_list,AFR_list):
		fields = line.strip().split("\t")
		self.seg_len += int(fields[2]) - int(fields[1])
		self.sample_num += 1
		self.temp_line_list.append(line)

		sample_name = fields[4]
		if sample_name in AFR_list:
			self.afr_num += 1
		if sample_name in LWK_list:
			self.lwk_num += 1
			self.lwk_len += int(fields[6]) - int(fields[5])

	def check_if_ILS(self):
		afr_ratio = self.afr_num / self.sample_num
		lwk_len_ratio = self.lwk_len / self.seg_len
		if afr_ratio >= 0.6 or ( self.lwk_num >= 3 and lwk_len_ratio >= 0.5):
			# print("ILS tract detected: AFR ratio = {}, LWK len ratio = {}".format(afr_ratio, lwk_len_ratio))
			return True
		else:
			return False

def get_sample_list(file_path):
	sample_list = []
	with open(file_path,"r") as f:
		for line in f:
			line = line.strip()
			if line and not line.startswith("#"):
				sample_list.append(line)
	return sample_list


if __name__ == '__main__':
	if len(sys.argv) != 4:
		print("Warning: Incorrect number of arguments.")
		print(f"Usage: python {sys.argv[0]} <total_bed> <total_detail> <AFR_sample_list> <LWK_sample_list> <out_prefix>")
		sys.exit(1)

	total_bed = sys.argv[1]
	total_detail = sys.argv[2]
	AFR_file = sys.argv[3]
	LWK_file = sys.argv[4]
	out_rmILS_prefix = sys.argv[5]

	'''
		Store AFR and LWK info
	'''
	LWK_list = get_sample_list(LWK_file)
	AFR_list = get_sample_list(AFR_file)


	'''
		Read Detail file, Determine wether the tract is ILS or not
	'''
	tract_dict = dict() 
	# key: <chrom> ; value: {'pos':ILS_flag}
	
	out_rmILS_detail = out_rmILS_prefix + ".detail"
	start_flag = 1
	ILS_flag_seg_num = 0
	with rich.progress.open(total_detail,"r",description="Process Detail File...") as detail_f, open(out_rmILS_detail,"w") as out_f:
		for line in detail_f:
			if line.startswith("#"):
				out_f.write(line)
				continue
			if start_flag:
				pre_chr = line.strip().split("\t")[0]
				pre_start_pos = int(line.strip().split("\t")[1])
				pre_tract = MergedTract()
				start_flag = 0
			cur_chr = line.strip().split("\t")[0]
			cur_start_pos = int(line.strip().split("\t")[1])
			if cur_chr == pre_chr and cur_start_pos == pre_start_pos:
				pre_tract.add_line(line,LWK_list,AFR_list)
			else:
				ILS_flag = pre_tract.check_if_ILS()
				if pre_chr not in tract_dict:
					tract_dict[pre_chr] = dict()
				tract_dict[pre_chr][pre_start_pos] = ILS_flag
				if not ILS_flag:
					out_f.write("".join(pre_tract.temp_line_list))
				else:
					ILS_flag_seg_num += 1
					#print("Non-ILS tract detected: {} {}".format(pre_chr, pre_start_pos))
				pre_chr = cur_chr
				pre_start_pos = cur_start_pos
				pre_tract = MergedTract()
				pre_tract.add_line(line,LWK_list,AFR_list)
		
		# Handle the last tract
		ILS_flag = pre_tract.check_if_ILS()
		if pre_chr not in tract_dict:
			tract_dict[pre_chr] = dict()
		tract_dict[pre_chr][pre_start_pos] = ILS_flag
		if not ILS_flag:
			out_f.write("".join(pre_tract.temp_line_list))
		else:
			ILS_flag_seg_num += 1
		
		print("Total ILS segments detected: {}".format(ILS_flag_seg_num))

	out_rmILS_bed = out_rmILS_prefix + ".bed"
	out_rmILS_bed_dict = dict()
	out_ILS_list = out_rmILS_prefix + ".ILS.list"
	with rich.progress.open(total_bed,"r",description="Remove the ILS in Bed...") as input_f, \
		open(out_rmILS_bed,"w") as out_bed_f, open(out_ILS_list,"w") as out_ILS_f, \
		open(out_rmILS_prefix + ".pickle","wb") as out_pickle_f:
		for line in input_f:
			if line.startswith("#"):
				out_bed_f.write(line)
				out_ILS_f.write(line)
				continue
			fields = line.strip().split("\t")
			chrom = fields[0]
			start_pos = int(fields[1])
			end_pos = int(fields[2])
			if chrom in tract_dict and start_pos in tract_dict[chrom]:
				if tract_dict[chrom][start_pos]: # ILS
					out_ILS_f.write(line)
				else: # Not-ILS
					out_bed_f.write(line)
					if chrom not in out_rmILS_bed_dict:
						out_rmILS_bed_dict[chrom] = []
					out_rmILS_bed_dict[chrom].append((start_pos, end_pos))
			else:
				print("Warning: {} {} not found in detail file".format(chrom, start_pos))
		pickle.dump(out_rmILS_bed_dict, out_pickle_f)

				

				

