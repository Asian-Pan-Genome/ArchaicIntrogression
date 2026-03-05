import sys
import pandas as pd
import pysam

df = pd.read_csv(sys.argv[1], sep="\t", header=0)
df = df[df['fisher_pvalue'] < 5e-8]
df = df[(df['REF'].str.len() == 1) & (df['ALT'].str.len() == 1)]
df['intro_allele'] = df.apply(lambda row: 'REF' if row['introgressed_0'] / (row['introgressed_0'] + row['introgressed_1']) > row['non_introgressed_0'] / (row['non_introgressed_0'] + row['non_introgressed_1']) else 'ALT', axis=1)
df['status'] = 'Low-confidence'

vcf_in = pysam.VariantFile(sys.argv[2])
for index, row in df.iterrows():
    chrom = row['chrom']
    pos = row['pos']
    ref = row['REF'].upper()
    alt = row['ALT'].upper()
    records = list(vcf_in.fetch(chrom, pos-1, pos))
    if row['intro_allele'] == 'REF':
        if len(records) == 0:
            df.at[index, 'status'] = 'High-confidence'
        else:
            for record in records:
                if record.pos == pos and record.ref.upper() == ref and alt == record.alts[0].upper():
                    if df.at[index, 'status'] == 'Low-confidence':
                        for sample, value_dict in record.samples.items():
                            if value_dict['GT'] in [(0,0), (0,1), (1,0)]:
                                df.at[index, 'status'] = 'High-confidence'
                                break
    else:
        for record in records:
            if record.pos == pos and record.ref.upper() == ref and alt == record.alts[0].upper():
                df.at[index, 'status'] = 'High-confidence'
                break
df = df[['chrom', 'pos', 'block_id', 'fisher_pvalue', 'REF', 'ALT', 'status']]
df.to_csv(sys.stdout, sep="\t", index=False, header=True)
vcf_in.close()
