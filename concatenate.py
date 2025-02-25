import pandas as pd

localDir = null
gbm_file = localDir + '/Gliobalstoma-Multiforme_tcga-rnaseq_counts.tsv'
lgg_file = localDir + '/Glioma_TCGA_counts.tcga-rnaseq.tsv'
output_file = '/concatenate.tsv'
#combine all tsv files in the list
combined_file = pd.concat([pd.read_csv(f) for f in [gbm_file, lgg_file]], sep='\t')

#export to tsv
combined_file.to_csv( output_file, sep='\t', index=False, encoding='utf-8-sig')