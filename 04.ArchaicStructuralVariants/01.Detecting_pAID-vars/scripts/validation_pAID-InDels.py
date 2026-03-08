import sys
import pandas as pd
import subprocess
import pysam

if len(sys.argv) != 7:
    print(f"Usage: python3 {sys.argv[0]} <input.freq_test.tsv> <input_A_ins.vcf> <input_A_del.vcf> <input_B_ins.vcf> <input_B_del.vcf> <output.tsv>")
    sys.exit(1)


df = pd.read_csv(sys.argv[1], sep="\t", header=0)
df = df.dropna(subset=['fisher_pvalue'])
df = df[df['fisher_pvalue'] < 5e-8]
#df['archaic_is_ref'] = df['introgressed_0'] > df['introgressed_1']
## actually, the result should be `introgressed_0` / (`introgressed_0` + `introgressed_1`) > `non_introgressed_0` / (`non_introgressed_0` + `non_introgressed_1`)
df['archaic_is_ref'] = (df['introgressed_0'] / (df['introgressed_0'] + df['introgressed_1'])) > (df['non_introgressed_0'] / (df['non_introgressed_0'] + df['non_introgressed_1']))


df_A = df[~df['archaic_is_ref']].copy()
df_B = df[df['archaic_is_ref']].copy()
df_A_ins = df_A[((df_A['REF'].str.len() == 1) & (df_A['ALT'].str.len() > 1)) & ((pd.isna(df_A['SVLEN'])) | ((df_A['SVLEN'] > -50) & (df_A['SVLEN'] < 50)))].copy()
df_A_del = df_A[((df_A['ALT'].str.len() == 1) & (df_A['REF'].str.len() > 1)) & ((pd.isna(df_A['SVLEN'])) | ((df_A['SVLEN'] > -50) & (df_A['SVLEN'] < 50)))].copy()
df_B_ins = df_B[((df_B['REF'].str.len() == 1) & (df_B['ALT'].str.len() > 1)) & ((pd.isna(df_B['SVLEN'])) | ((df_B['SVLEN'] > -50) & (df_B['SVLEN'] < 50)))].copy()
df_B_del = df_B[((df_B['ALT'].str.len() == 1) & (df_B['REF'].str.len() > 1)) & ((pd.isna(df_B['SVLEN'])) | ((df_B['SVLEN'] > -50) & (df_B['SVLEN'] < 50)))].copy()


df_A_ins['status'] = 'Low-confidence'
df_A_del['status'] = 'Low-confidence'
with pysam.VariantFile(sys.argv[2]) as vcf_in:
    for rec in vcf_in.fetch():
        _, value_dict = list(rec.samples.items())[0]
        if value_dict['BD'] == 'TP':
            chrom = rec.chrom
            pos = rec.pos
            ref = rec.ref.upper()
            alt = rec.alts[0].upper()
            df_A_ins.loc[(df_A_ins['chrom'] == chrom) & (df_A_ins['pos'] == pos) & (df_A_ins['REF'].str.upper() == ref) & (df_A_ins['ALT'].str.upper() == alt), 'status'] = 'High-confidence'
df_A_ins = df_A_ins[['chrom', 'pos', 'block_id', 'fisher_pvalue', 'REF', 'ALT', 'status']]

with pysam.VariantFile(sys.argv[3]) as vcf_in:
    for rec in vcf_in.fetch():
        _, value_dict = list(rec.samples.items())[0]
        if value_dict['BD'] == 'TP':
            chrom = rec.chrom
            pos = rec.pos
            ref = rec.ref.upper()
            alt = rec.alts[0].upper()
            df_A_del.loc[(df_A_del['chrom'] == chrom) & (df_A_del['pos'] == pos) & (df_A_del['REF'].str.upper() == ref) & (df_A_del['ALT'].str.upper() == alt), 'status'] = 'High-confidence'
df_A_del = df_A_del[['chrom', 'pos', 'block_id', 'fisher_pvalue', 'REF', 'ALT', 'status']]
df_A_final = pd.concat([df_A_ins, df_A_del], axis=0)
#df_A_final.to_csv(sys.argv[4], sep="\t", index=False, header=True)


df_B_ins['status'] = 'High-confidence'
df_B_del['status'] = 'High-confidence'
with pysam.VariantFile(sys.argv[4]) as vcf_in:
    for rec in vcf_in.fetch():
        _, value_dict = list(rec.samples.items())[0]
        if value_dict['BD'] == 'TP':
            chrom = rec.chrom
            pos = rec.pos
            ref = rec.ref.upper()
            alt = rec.alts[0].upper()
            df_B_ins.loc[(df_B_ins['chrom'] == chrom) & (df_B_ins['pos'] == pos) & (df_B_ins['REF'].str.upper() == ref) & (df_B_ins['ALT'].str.upper() == alt), 'status'] = 'Low-confidence'
df_B_ins = df_B_ins[['chrom', 'pos', 'block_id', 'fisher_pvalue', 'REF', 'ALT', 'status']]

with pysam.VariantFile(sys.argv[5]) as vcf_in:
    for rec in vcf_in.fetch():
        _, value_dict = list(rec.samples.items())[0]
        if value_dict['BD'] == 'TP':
            chrom = rec.chrom
            pos = rec.pos
            ref = rec.ref.upper()
            alt = rec.alts[0].upper()
            df_B_del.loc[(df_B_del['chrom'] == chrom) & (df_B_del['pos'] == pos) & (df_B_del['REF'].str.upper() == ref) & (df_B_del['ALT'].str.upper() == alt), 'status'] = 'Low-confidence'
df_B_del = df_B_del[['chrom', 'pos', 'block_id', 'fisher_pvalue', 'REF', 'ALT', 'status']]
df_B_final = pd.concat([df_B_ins, df_B_del], axis=0)

df_final = pd.concat([df_A_final, df_B_final], axis=0)
df_final.to_csv(sys.argv[6], sep="\t", index=False, header=True)
