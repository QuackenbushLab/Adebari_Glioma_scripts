import pandas as pd
#Pandas is a packages for working with datasets
#define file path first 
gbm_file_path = null
lgg_file_path = null
split_file_dir = null
#Read files into dataframes
gbm_df = pd.read_csv(gbm_file_path, sep='\t')
lgg_df = pd.read_csv(lgg_file_path, sep='\t')
#Filter dataframes
gbm_male = gbm_df[gbm_df['gender'] == 'male']
gbm_female = gbm_df[gbm_df['gender'] == 'female']
lgg_male = lgg_df[lgg_df['gender'] == 'male']
lgg_female = lgg_df[lgg_df['gender'] == 'female']
#Overlapping genes
#Path to new files
gbm_male_path = split_file_dir + '/GBM_Male.tsv'
gbm_female_path = split_file_dir + '/GBM_Female.tsv'
lgg_male_path = split_file_dir + '/LGG_Male.tsv'
lgg_female_path = split_file_dir + '/LGG_Female.tsv'
#save files
gbm_male.to_csv(gbm_male_path, sep='\t', index=False)
gbm_female.to_csv(gbm_female_path, sep='\t', index=False)
lgg_male.to_csv(lgg_male_path, sep='\t', index=False)
lgg_female.to_csv(lgg_female_path, sep='\t', index=False)
print("Files have been sucessfully seperated and saved on your desktop!")