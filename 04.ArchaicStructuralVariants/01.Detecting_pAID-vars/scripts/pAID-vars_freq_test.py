import sys
import pandas as pd
import pysam
import scipy.stats as stats
import numpy as np
import re

if len(sys.argv) != 4:
    print(f"Usage: python {sys.argv[0]} input.detail gene input.INV.tsv")
    sys.exit(1)


input_detail = pd.read_csv(sys.argv[1], sep='\t', header=None, names=['chr', 'start', 'end', 'block_id', 'sample'])
gene = sys.argv[2]
INV_tsv = pd.read_csv(sys.argv[3], sep='\t', header=0)

print('chrom\tpos\tblock_id\tchi_square_pvalue\tfisher_pvalue\tintrogressed_0\tintrogressed_1\tintrogressed_missing\tnon_introgressed_0\tnon_introgressed_1\tnon_introgressed_missing\tSVTYPE\tSVLEN\tREF\tALT')
with pysam.VariantFile(f'{gene}.vcf.gz') as vcf:
    for block_id, group in input_detail.groupby('block_id'):
        if len(group['sample']) < 10:
            continue
        chrom = group['chr'].unique()[0]
        start = group['start'].unique()[0]
        end = group['end'].unique()[0]


        for rec in vcf.fetch(contig=chrom, start=start, stop=end):
            introgressed = {'0':0, '1':0, 'missing':0}
            non_introgressed = {'0':0, '1':0, 'missing':0}
            
            for sample, value_dict in rec.samples.items():
                if len(value_dict['GT']) == 2:
                    if f'{sample}_hap1' in group['sample'].values:
                        if value_dict['GT'][0] is None:
                            introgressed['missing'] += 1
                        elif value_dict['GT'][0] == 0:
                            introgressed['0'] += 1
                        elif value_dict['GT'][0] == 1:
                            introgressed['1'] += 1
                    else:
                        if value_dict['GT'][0] is None:
                            non_introgressed['missing'] += 1
                        elif value_dict['GT'][0] == 0:
                            non_introgressed['0'] += 1
                        elif value_dict['GT'][0] == 1:
                            non_introgressed['1'] += 1
                    
                    if f'{sample}_hap2' in group['sample'].values:
                        if value_dict['GT'][1] is None:
                            introgressed['missing'] += 1
                        elif value_dict['GT'][1] == 0:
                            introgressed['0'] += 1
                        elif value_dict['GT'][1] == 1:
                            introgressed['1'] += 1
                    else:
                        if value_dict['GT'][1] is None:
                            non_introgressed['missing'] += 1
                        elif value_dict['GT'][1] == 0:
                            non_introgressed['0'] += 1
                        elif value_dict['GT'][1] == 1:
                            non_introgressed['1'] += 1
            

            if introgressed['0'] + non_introgressed['0'] == 0 or introgressed['1'] + non_introgressed['1'] == 0:
                chi_p = np.nan
                fisher_p = np.nan
            else:
                if introgressed['missing'] == 0 and non_introgressed['missing'] == 0:
                    chi_p = stats.chi2_contingency([[introgressed['0'], introgressed['1']], [non_introgressed['0'], non_introgressed['1']]]).pvalue
                else:
                    try:
                        chi_p = stats.chi2_contingency([[introgressed['0'], introgressed['1'], introgressed['missing']], [non_introgressed['0'], non_introgressed['1'], non_introgressed['missing']]]).pvalue
                    except:
                        sys.exit(f'{rec.chrom}\t{rec.pos}\t{introgressed["0"]}\t{introgressed["1"]}\t{introgressed["missing"]}\t{non_introgressed["0"]}\t{non_introgressed["1"]}\t{non_introgressed["missing"]}\t{rec.info.get("SVTYPE", "")}\t{rec.info.get("SVLEN", ("", ""))[0]}\t{rec.alleles[0]}\t{rec.alleles[1]}')
            
                fisher_p = stats.fisher_exact([[introgressed['0'], introgressed['1']], [non_introgressed['0'], non_introgressed['1']]]).pvalue
            
            print(f'{rec.chrom}\t{rec.pos}\t{block_id}\t{chi_p}\t{fisher_p}\t{introgressed["0"]}\t{introgressed["1"]}\t{introgressed["missing"]}\t{non_introgressed["0"]}\t{non_introgressed["1"]}\t{non_introgressed["missing"]}\t{rec.info.get("SVTYPE", "")}\t{rec.info.get("SVLEN", ("", ""))[0]}\t{rec.alleles[0]}\t{rec.alleles[1]}')

        
        tmp_INV_tsv = INV_tsv[(INV_tsv['CHR'] == chrom) & (INV_tsv['START'] >= start) & (INV_tsv['END'] <= end)]
        if not tmp_INV_tsv.empty:
            for _, row in tmp_INV_tsv.iterrows():
                introgressed = {'0':0, '1':0, 'missing':0}
                non_introgressed = {'0':0, '1':0, 'missing':0}
                for col in tmp_INV_tsv.columns[5:]:
                    gt = str(row[col])
                    col = col.replace('.mat', '_hap1').replace('.pat', '_hap2')
                    if col in group['sample'].values:
                        if gt is None or gt == '.':
                            introgressed['missing'] += 1
                        elif gt == '0':
                            introgressed['0'] += 1
                        elif gt == '1':
                            introgressed['1'] += 1
                    else:
                        if gt is None or gt == '.':
                            non_introgressed['missing'] += 1
                        elif gt == '0':
                            non_introgressed['0'] += 1
                        elif gt == '1':
                            non_introgressed['1'] += 1

                if introgressed['0'] + non_introgressed['0'] == 0 or introgressed['1'] + non_introgressed['1'] == 0:
                    chi_p = np.nan
                    fisher_p = np.nan
                else:
                    if introgressed['missing'] == 0 and non_introgressed['missing'] == 0:
                        chi_p = stats.chi2_contingency([[introgressed['0'], introgressed['1']], [non_introgressed['0'], non_introgressed['1']]]).pvalue
                    else:
                        chi_p = stats.chi2_contingency([[introgressed['0'], introgressed['1'], introgressed['missing']], [non_introgressed['0'], non_introgressed['1'], non_introgressed['missing']]]).pvalue
                    
                    fisher_p = stats.fisher_exact([[introgressed['0'], introgressed['1']], [non_introgressed['0'], non_introgressed['1']]]).pvalue
                print(f'{row["CHR"]}\t{row["START"]}\t{block_id}\t{chi_p}\t{fisher_p}\t{introgressed["0"]}\t{introgressed["1"]}\t{introgressed["missing"]}\t{non_introgressed["0"]}\t{non_introgressed["1"]}\t{non_introgressed["missing"]}\tINV\t{row["LENGTH"]}\t""\t""')
