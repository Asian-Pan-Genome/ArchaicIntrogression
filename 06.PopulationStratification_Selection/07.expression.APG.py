import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
from matplotlib.ticker import AutoMinorLocator
import matplotlib.patches as mpatches
from scipy.stats import linregress
import numpy as np
from statsmodels.regression.linear_model import OLS
from statsmodels.stats.multitest import multipletests
import pysam
import sys
import subprocess
from matplotlib.backends.backend_pdf import PdfPages # Import multipage PDF backend

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42

if len(sys.argv) != 6:
    sys.exit(f'Usage: python {sys.argv[0]} <input_sv_tsv> <input_TMM> <gene_id> <input.vcf> <output_prefix>')

# 1. Pre-read expression data (read only once)
print("Loading Expression Data...")
expression_df_raw = pd.read_csv(sys.argv[2], sep='\t', header=0)
expression_df_raw = expression_df_raw.set_index('Gene')
platform = sys.argv[2].split('.')[0].split('_')[-1]
gene_id = sys.argv[3]

if platform == 'ONT':
    gene_id = gene_id.split('.')[0]

# Extract expression data for specific gene as baseline data
if gene_id not in expression_df_raw.index:
    sys.exit(f'Error: Gene ID {gene_id} not found in expression matrix.')

base_df = expression_df_raw.loc[gene_id].to_frame('TMM').iloc[3:, :]
# Ensure TMM column is numeric
base_df['TMM'] = pd.to_numeric(base_df['TMM'])

# 2. Index VCF (execute only once)
print("Indexing VCF...")
subprocess.run('bcftools index -f {}'.format(sys.argv[4]), shell=True, check=True)

input_sv_file = sys.argv[1]
input_vcf_file = sys.argv[4]
output_prefix = sys.argv[5]

# 3. Open PDF and VCF files, start loop processing
print(f"Processing sites from {input_sv_file}...")

with PdfPages(f'{output_prefix}.pdf') as pdf, \
     pysam.VariantFile(input_vcf_file) as vcf_in, \
     open(input_sv_file, 'r') as f:

    for line_idx, line in enumerate(f):
        line = line.strip()
        if not line: continue
        
        try:
            # Parse data from each line
            tmp = line.split('\t')
            chrom = tmp[0]
            pos = int(tmp[1])
            introgressed_0 = int(tmp[5])
            introgressed_1 = int(tmp[6])
            non_introgressed_0 = int(tmp[8])
            non_introgressed_1 = int(tmp[9])
            ref_allele = tmp[13]
            alt_allele = tmp[14]

            # Determine intro_allele
            if introgressed_0 / (introgressed_0 + introgressed_1) > non_introgressed_0 / (non_introgressed_0 + non_introgressed_1):
                intro_allele = ref_allele
            else:
                intro_allele = alt_allele

            # Create a new DataFrame copy for current site
            curr_df = base_df.copy()
            curr_df['Genotype'] = np.nan

            # Find variant in VCF
            all_recs = []
            # Note: pysam fetch uses 0-based coordinates, and end is exclusive
            try:
                records = vcf_in.fetch(chrom, pos-1, pos-1+len(ref_allele))
            except ValueError as e:
                print(f"Warning: Could not fetch region {chrom}:{pos}. Error: {e}")
                continue

            for rec in records:
                if rec.pos == pos:
                    # Match ref/alt
                    if rec.ref.upper() == ref_allele.upper() and rec.alts[0].upper() == alt_allele.upper():
                        if intro_allele == ref_allele:
                            all_recs.append((rec, 0)) # intro is ref (0)
                        else:
                            all_recs.append((rec, 1)) # intro is alt (1)
                    elif rec.ref.upper() == alt_allele.upper() and rec.alts[0].upper() == ref_allele.upper():
                        if intro_allele == ref_allele:
                            all_recs.append((rec, 1)) 
                        else:
                            all_recs.append((rec, 0)) 

            if len(all_recs) == 0:
                print(f'Warning: No matching variant found in VCF for {chrom}:{pos}. Skipping.')
                continue
            elif len(all_recs) > 1:
                print(f'Warning: Multiple matching variants found for {chrom}:{pos}. Skipping.')
                continue
            
            # Extract genotypes
            rec, intro_allele_index = all_recs[0]
            
            for sample, value_dict in rec.samples.items():
                if sample in curr_df.index:
                    genotype = value_dict['GT']
                    if genotype is None or None in genotype:
                        continue # Skip samples with missing GT
                    
                    if genotype[0] == intro_allele_index and genotype[1] == intro_allele_index:
                        curr_df.at[sample, 'Genotype'] = 'intro/intro'
                    elif genotype[0] != intro_allele_index and genotype[1] != intro_allele_index:
                        curr_df.at[sample, 'Genotype'] = 'non-intro/non-intro'
                    else:
                        curr_df.at[sample, 'Genotype'] = 'intro/non-intro'

            # Clean data again, remove samples without matching genotypes
            curr_df = curr_df.dropna(subset=['Genotype'])
            
            # Skip statistics if too few samples
            if len(curr_df) < 3:
                print(f"Warning: Not enough samples with genotype info for {chrom}:{pos}. Skipping.")
                continue

            curr_df['Genotype_numeric'] = curr_df['Genotype'].astype('category').cat.codes

            # Check if there is sufficient variation for regression
            if len(curr_df['Genotype'].unique()) < 2:
                print(f"Warning: Only one genotype class found for {chrom}:{pos}. Skipping OLS.")
                p_value_text = "P = N/A (Single Genotype)"
            else:
                # Statistical analysis (OLS)
                results = OLS(curr_df['TMM'], curr_df['Genotype_numeric']).fit()
                p_value_text = f'P = {results.pvalues.iloc[0]:.2e}'
            
            # Plotting
            color_map = plt.cm.get_cmap('Set1')
            fig, ax = plt.subplots(figsize=(10, 8))
            
            # Define order
            order_list = ['non-intro/non-intro', 'intro/non-intro', 'intro/intro']
            
            sns.boxplot(x='Genotype', y='TMM', data=curr_df, 
                        order=order_list, 
                        hue='Genotype', 
                        palette=[color_map(0), color_map(1), color_map(2)], 
                        showcaps=False, showfliers=False, dodge=False, ax=ax)
            
            sns.swarmplot(x='Genotype', y='TMM', data=curr_df, 
                          order=order_list, 
                          color='lightgray', alpha=0.7, ax=ax)

            # Add P-value
            ax.text(0.05, 0.95, p_value_text, transform=ax.transAxes, fontsize=14, verticalalignment='top')
            
            # Calculate sample counts and update x-axis labels
            counts = [len(curr_df[curr_df['Genotype'] == cat]) for cat in order_list]
            new_xticklabels = [f'{cat}\n($n$={n})' for cat, n in zip(order_list, counts)]
            ax.set_xticklabels(new_xticklabels)

            ax.set_xlabel('Allele', fontsize=16)
            ax.set_ylabel('TMM', fontsize=16)
            # Add site information in title for differentiation
            ax.set_title(f'APGp1 ({platform}) - {chrom}:{pos}', fontsize=18)
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            sns.despine()
            plt.tight_layout()
            
            # Save current figure to PDF
            pdf.savefig(fig)
            plt.close(fig) # Close figure to free memory
            
            print(f"Processed {chrom}:{pos}")

        except Exception as e:
            print(f"Error processing line {line_idx+1} ({line.split()[0]}:{line.split()[1] if len(line.split())>1 else '?'}): {e}")
            plt.close('all') # Ensure to close plots when error occurs
            continue

print(f"Done. Output saved to {output_prefix}.pdf")