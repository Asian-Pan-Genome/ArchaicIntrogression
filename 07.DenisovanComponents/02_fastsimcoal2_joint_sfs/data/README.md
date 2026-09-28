# Data included in this directory

`observed_sfs/YRI10_CHB20_DEN_1SNPper100kb_MSFS.obs` is the folded joint SFS used for the eight-model fastsimcoal2 comparison. It contains 27,160 polymorphic sites and chromosome sample sizes of 20 YRI, 40 CHB and 2 Denisova 3 chromosomes.

The autosomes were divided into non-overlapping 100-kb blocks, and one SNP was selected from each non-empty block using reservoir sampling with a fixed seed. Folding was based on the allele with the lower total count across all 62 chromosomes. When the total allele count was exactly 31, the site count was divided equally between the two complementary entries; this accounts for the fractional values in a small number of SFS cells.

The `samples/` directory records the 10 YRI and 20 CHB individuals used to construct the SFS. Raw modern-human and Denisovan VCF files are not redistributed here.
