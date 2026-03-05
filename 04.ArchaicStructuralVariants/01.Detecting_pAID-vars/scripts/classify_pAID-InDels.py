import sys
import pandas as pd
import subprocess
import pysam

if len(sys.argv) != 6:
    print(f"Usage: python3 {sys.argv[0]} <input.freq_test.tsv> <output_A_ins_vcf> <output_A_del_vcf> <output_B_ins_vcf> <output_B_del_vcf>")
    sys.exit(1)


df = pd.read_csv(sys.argv[1], sep="\t", header=0)
df = df.dropna(subset=['fisher_pvalue'])
df = df[df['fisher_pvalue'] < 5e-8]
df['archaic_is_ref'] = (df['introgressed_0'] / (df['introgressed_0'] + df['introgressed_1'])) > (df['non_introgressed_0'] / (df['non_introgressed_0'] + df['non_introgressed_1']))

df_A = df[~df['archaic_is_ref']].copy()
df_B = df[df['archaic_is_ref']].copy()
df_A_ins = df_A[((df_A['REF'].str.len() == 1) & (df_A['ALT'].str.len() > 1)) & ((pd.isna(df_A['SVLEN'])) | ((df_A['SVLEN'] > -50) & (df_A['SVLEN'] < 50)))].copy()
df_A_del = df_A[((df_A['ALT'].str.len() == 1) & (df_A['REF'].str.len() > 1)) & ((pd.isna(df_A['SVLEN'])) | ((df_A['SVLEN'] > -50) & (df_A['SVLEN'] < 50)))].copy()
df_B_ins = df_B[((df_B['REF'].str.len() == 1) & (df_B['ALT'].str.len() > 1)) & ((pd.isna(df_B['SVLEN'])) | ((df_B['SVLEN'] > -50) & (df_B['SVLEN'] < 50)))].copy()
df_B_del = df_B[((df_B['ALT'].str.len() == 1) & (df_B['REF'].str.len() > 1)) & ((pd.isna(df_B['SVLEN'])) | ((df_B['SVLEN'] > -50) & (df_B['SVLEN'] < 50)))].copy()



header_lines = [
    '##fileformat=VCFv4.2',
    '##FILTER=<ID=PASS,Description="All filters passed">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##contig=<ID=chr1,length=248387328>',
    '##contig=<ID=chr2,length=242696752>',
    '##contig=<ID=chr3,length=201105948>',
    '##contig=<ID=chr4,length=193574945>',
    '##contig=<ID=chr5,length=182045439>',
    '##contig=<ID=chr6,length=172126628>',
    '##contig=<ID=chr7,length=160567428>',
    '##contig=<ID=chr8,length=146259331>',
    '##contig=<ID=chr9,length=150617247>',
    '##contig=<ID=chr10,length=134758134>',
    '##contig=<ID=chr11,length=135127769>',
    '##contig=<ID=chr12,length=133324548>',
    '##contig=<ID=chr13,length=113566686>',
    '##contig=<ID=chr14,length=101161492>',
    '##contig=<ID=chr15,length=99753195>',
    '##contig=<ID=chr16,length=96330374>',
    '##contig=<ID=chr17,length=84276897>',
    '##contig=<ID=chr18,length=80542538>',
    '##contig=<ID=chr19,length=61707364>',
    '##contig=<ID=chr20,length=66210255>',
    '##contig=<ID=chr21,length=45090682>',
    '##contig=<ID=chr22,length=51324926>',
    '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'])
]
ins_lines = header_lines[:]
del_lines = header_lines[:]

for rowindex, row in df_A_ins.iterrows():
    chrom = row['chrom']
    pos = row['pos']
    ref = row['REF']
    alt = row['ALT']
    vcf_line = '\t'.join([str(chrom), str(pos), '.', ref, alt, '.', 'PASS', '.', 'GT', '1/1'])
    ins_lines.append(vcf_line)
with open(sys.argv[2], 'w') as f_out:
    f_out.write('\n'.join(ins_lines) + '\n')
subprocess.run(f'bcftools sort {sys.argv[2]} -Oz -o {sys.argv[2]}.gz', shell=True)
subprocess.run(f'rm -rf {sys.argv[2]}', shell=True)
subprocess.run(f'bcftools index {sys.argv[2]}.gz', shell=True)

for rowindex, row in df_A_del.iterrows():
    chrom = row['chrom']
    pos = row['pos']
    ref = row['REF']
    alt = row['ALT']
    vcf_line = '\t'.join([str(chrom), str(pos), '.', ref, alt, '.', 'PASS', '.', 'GT', '1/1'])
    del_lines.append(vcf_line)
with open(sys.argv[3], 'w') as f_out:
    f_out.write('\n'.join(del_lines) + '\n')
subprocess.run(f'bcftools sort {sys.argv[3]} -Oz -o {sys.argv[3]}.gz', shell=True)
subprocess.run(f'rm -rf {sys.argv[3]}', shell=True)
subprocess.run(f'bcftools index {sys.argv[3]}.gz', shell=True)

ins_lines = header_lines[:]
del_lines = header_lines[:]
for rowindex, row in df_B_ins.iterrows():
    chrom = row['chrom']
    pos = row['pos']
    ref = row['REF']
    alt = row['ALT']
    vcf_line = '\t'.join([str(chrom), str(pos), '.', ref, alt, '.', 'PASS', '.', 'GT', '1/1'])
    ins_lines.append(vcf_line)
with open(sys.argv[4], 'w') as f_out:
    f_out.write('\n'.join(ins_lines) + '\n')
subprocess.run(f'bcftools sort {sys.argv[4]} -Oz -o {sys.argv[4]}.gz', shell=True)
subprocess.run(f'rm -rf {sys.argv[4]}', shell=True)
subprocess.run(f'bcftools index {sys.argv[4]}.gz', shell=True)

for rowindex, row in df_B_del.iterrows():
    chrom = row['chrom']
    pos = row['pos']
    ref = row['REF']
    alt = row['ALT']
    vcf_line = '\t'.join([str(chrom), str(pos), '.', ref, alt, '.', 'PASS', '.', 'GT', '1/1'])
    del_lines.append(vcf_line)
with open(sys.argv[5], 'w') as f_out:
    f_out.write('\n'.join(del_lines) + '\n')
subprocess.run(f'bcftools sort {sys.argv[5]} -Oz -o {sys.argv[5]}.gz', shell=True)
subprocess.run(f'rm -rf {sys.argv[5]}', shell=True)
subprocess.run(f'bcftools index {sys.argv[5]}.gz', shell=True)

