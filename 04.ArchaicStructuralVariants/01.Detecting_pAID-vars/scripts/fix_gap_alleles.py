import sys
import pysam
from collections import OrderedDict

def filter_gap_alleles(vcf_in_path, vcf_out_path):
    """
    Filters alleles containing 'N' or 'n' from a VCF file and updates
    genotypes and allele-dependent INFO fields accordingly.
    """
    try:
        # 使用pysam打开输入和输出VCF文件
        with pysam.VariantFile(vcf_in_path) as f_in, \
             pysam.VariantFile(vcf_out_path, 'w', header=f_in.header) as f_out:

            for rec in f_in.fetch():
                alleles = rec.alleles
                indices_to_remove = set()

                # 1. 识别含有'N/n'的等位基因（跳过REF等位基因）
                for i, allele in enumerate(alleles):
                    if i > 0 and ('N' in allele or 'n' in allele):
                        indices_to_remove.add(i)

                # 如果没有需要移除的等位基因，直接写入并继续
                if not indices_to_remove:
                    f_out.write(rec)
                    continue
                if len(indices_to_remove) == len(alleles) - 1:
                    continue

                # 2. 创建旧索引到新索引的映射
                new_alleles = []
                index_map = []
                new_idx_counter = 0
                for i, allele in enumerate(alleles):
                    if i in indices_to_remove:
                        index_map.append(None) # 标记为删除
                    else:
                        new_alleles.append(allele)
                        index_map.append(new_idx_counter)
                        new_idx_counter += 1
                
                # 创建一个新的VCF记录
                new_info = dict(rec.info)
                new_info['AC'] = [ac for i, ac in enumerate(rec.info['AC'], start=1) if i not in indices_to_remove]
                new_info['AF'] = [af for i, af in enumerate(rec.info['AF'], start=1) if i not in indices_to_remove]
                new_info['AT'] = [at for i, at in enumerate(rec.info['AT']) if i not in indices_to_remove]
                new_info['AN'] = rec.info['AN'] - sum([rec.info['AC'][i-1] for i in indices_to_remove])
                new_rec = f_out.header.new_record(contig=rec.chrom, start=rec.start, stop=rec.stop, alleles=new_alleles, id=rec.id, qual=rec.qual, filter=rec.filter, info=new_info)
                

                # 4. 遍历并更新每个样本的基因型
                for sample_name in rec.samples:
                    original_sample = rec.samples[sample_name]
                    new_sample = new_rec.samples[sample_name]
                    
                    original_gt = original_sample['GT']
                    new_gt = []

                    for allele_idx in original_gt:
                        if allele_idx is None: # 本来就是缺失基因型
                            new_gt.append(None)
                            continue
                        
                        mapped_idx = index_map[allele_idx]
                        new_gt.append(mapped_idx)
                    
                    new_sample['GT'] = tuple(new_gt)
                    new_sample.phased = original_sample.phased

                # 写入修改后的记录
                f_out.write(new_rec)

    except FileNotFoundError:
        print(f"Error: Input file not found at '{vcf_in_path}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <input.vcf> <output.vcf>")
        sys.exit(1)
    
    input_vcf = sys.argv[1]
    output_vcf = sys.argv[2]
    
    print(f"Filtering VCF '{input_vcf}' for gap alleles...")
    filter_gap_alleles(input_vcf, output_vcf)
    print(f"Processing complete. Output written to '{output_vcf}'.")