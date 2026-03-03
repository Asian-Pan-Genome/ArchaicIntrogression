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
from matplotlib.backends.backend_pdf import PdfPages  # New: for multipage PDF

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42

# Check parameters
if len(sys.argv) != 4:
    sys.exit(f'Usage: python {sys.argv[0]} <input_sv_tsv> <gene_id> <output_prefix>')

# 1. Load expression data outside the loop (optimization: read file only once)
print("Loading expression data...", file=sys.stderr)
df_mage_raw = pd.read_csv('./TMM.filtered.TSS.MAGE.v1.0.bed', sep='\t', header=0)
df_mage_raw = df_mage_raw.set_index('ID')

# Check if gene exists
if sys.argv[2] not in df_mage_raw.index:
    sys.exit(f"Error: Gene {sys.argv[2]} not found in expression file.")

# Extract base data for specific gene (since it's the same gene, this part is also universal)
df_gene_base = df_mage_raw.loc[sys.argv[2]]
df_gene_base = df_gene_base.to_frame('TMM').iloc[3:, :]
# Convert to numeric in advance to avoid repeated conversion in the loop
df_gene_base['TMM'] = pd.to_numeric(df_gene_base['TMM'])

# 2. Open VCF and PDF objects outside the loop
vcf_path = '1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.3202.bcf.gz'
output_pdf_name = f'{sys.argv[3]}.pdf'

print(f"Processing variants and saving to {output_pdf_name}...", file=sys.stderr)

with pysam.VariantFile(vcf_path) as vcf_in, PdfPages(output_pdf_name) as pdf:
    
    with open(sys.argv[1], 'r') as f:
        for line_idx, line in enumerate(f):
            # Remove empty lines
            if not line.strip():
                continue
                
            tmp = line.strip().split('\t')
            chrom = tmp[0]
            pos = int(tmp[1])
            introgressed_0 = int(tmp[5])
            introgressed_1 = int(tmp[6])
            non_introgressed_0 = int(tmp[8])
            non_introgressed_1 = int(tmp[9])
            SV_type = tmp[11]
            ref_allele = tmp[13]
            alt_allele = tmp[14]
            
            # Logic remains the same here
            if introgressed_0 / (introgressed_0 + introgressed_1) > non_introgressed_0 / (non_introgressed_0 + non_introgressed_1):
                intro_allele = ref_allele
            else:
                intro_allele = alt_allele

            # 3. Create a new DataFrame copy for the current site
            # Must copy, otherwise Genotype from the previous loop will contaminate the current loop
            df_mage = df_gene_base.copy()
            df_mage['Genotype'] = np.nan

            # VCF search logic
            all_recs = []
            true_exit = False
            
            # Note: pysam fetch interval is left-closed, right-open, and 0-based, keeping your original logic
            try:
                # Add try-except to prevent fetch from going out of bounds or chromosome non-existence from interrupting the program
                iter_vcf = vcf_in.fetch(chrom, pos-1, pos-1+len(ref_allele))
            except ValueError as e:
                print(f"Skipping {chrom}:{pos} - {e}", file=sys.stderr)
                continue

            for rec in iter_vcf:
                if rec.pos == pos:
                    if rec.ref.upper() == ref_allele.upper() and rec.alts[0].upper() == alt_allele.upper():
                        if intro_allele == ref_allele:
                            all_recs.append((rec, 0))
                        else:
                            all_recs.append((rec, 1))
                        true_exit = True
                        break
                    elif rec.ref.upper() == alt_allele.upper() and rec.alts[0].upper() == ref_allele.upper():
                        if intro_allele == ref_allele:
                            all_recs.append((rec, 1))
                        else:
                            all_recs.append((rec, 0))
                        true_exit = True
                        break
            
            # Logic modification note: Change sys.exit to continue to avoid single line error causing full exit
            if true_exit == False:
                print(f'Warning: no matching variant found in VCF for position {chrom}:{pos}. Skipping.', file=sys.stderr)
                continue
                    
            if len(all_recs) == 0:
                print(f'Warning: no matching variant found (empty recs) for position {chrom}:{pos}. Skipping.', file=sys.stderr)
                continue
            elif len(all_recs) > 1:
                print(f'Warning: multiple matching variants found for position {chrom}:{pos}. Skipping.', file=sys.stderr)
                continue
            else:
                rec, intro_allele_index = all_recs[0]
                for sample, value_dict in rec.samples.items():
                    if sample in df_mage.index:
                        genotype = value_dict['GT']
                        if None in genotype:
                            # Similarly, skip this sample or issue warning instead of exiting
                            print(f'Warning: unexpected genotype for sample {sample} at {chrom}:{pos}', file=sys.stderr)
                            continue
                        elif genotype[0] == intro_allele_index and genotype[1] == intro_allele_index:
                            df_mage.at[sample, 'Genotype'] = 'intro/intro'
                        elif genotype[0] != intro_allele_index and genotype[1] != intro_allele_index:
                            df_mage.at[sample, 'Genotype'] = 'non-intro/non-intro'
                        else:
                            df_mage.at[sample, 'Genotype'] = 'intro/non-intro'
            
            # Clean data
            df_mage = df_mage.dropna(subset=['Genotype'])
            
            # If no samples remain after filtering, skip plotting
            if df_mage.empty:
                print(f"Warning: No valid samples left for {chrom}:{pos}. Skipping.", file=sys.stderr)
                continue

            df_mage['Genotype_numeric'] = df_mage['Genotype'].astype('category').cat.codes
            
            # 4. Plotting logic
            color_map = plt.cm.get_cmap('Set1')
            fig, ax = plt.subplots(figsize=(10, 8))
            
            order_list = ['non-intro/non-intro', 'intro/non-intro', 'intro/intro']

            sns.boxplot(x='Genotype', y='TMM', data=df_mage, 
                        order=order_list, 
                        hue='Genotype', 
                        palette=[color_map(0), color_map(1), color_map(2)], 
                        showcaps=False, showfliers=False, dodge=False, ax=ax)
            
            sns.swarmplot(x='Genotype', y='TMM', data=df_mage, 
                          order=order_list, 
                          color='lightgray', alpha=0.7, ax=ax)

            # Statistical test
            try:
                results = OLS(df_mage['TMM'], df_mage['Genotype_numeric']).fit()
                p_val_text = f'P = {results.pvalues.iloc[0]:.2e}'
            except Exception as e:
                p_val_text = 'P = N/A' # Prevent regression failure (e.g., too little data)
            
            ax.text(0.05, 0.95, p_val_text, transform=ax.transAxes, fontsize=14, verticalalignment='top')
            
            # Calculate sample counts and update x-axis labels
            counts = [len(df_mage[df_mage['Genotype'] == cat]) for cat in order_list]
            new_xticklabels = [f'{cat}\n($n$={n})' for cat, n in zip(order_list, counts)]
            ax.set_xticklabels(new_xticklabels)

            ax.set_xlabel('Allele', fontsize=16)
            plt.ylabel('TMM', fontsize=16)
            
            # Modify title: add position information to distinguish in multipage PDF
            ax.set_title(f'MAGE - {chrom}:{pos}', fontsize=18)
            
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            sns.despine()
            plt.tight_layout()
            
            # 5. Save to PDF object instead of directly to file
            pdf.savefig(fig)
            plt.close(fig) # Close current figure to free memory
            
            print(f"Processed {chrom}:{pos}", file=sys.stderr)

print("Done.")