import pandas as pd
#file paths
split_file_dir = null
matrix_dir = null
LGG_Female_path =split_file_dir + '/LGG_Female.tsv'
LGG_Male_path =split_file_dir + '/LGG_Male.tsv'
GBM_Female_path =split_file_dir + '/GBM_Female.tsv'
GBM_Male_path = split_file_dir + '/GBM_Male.tsv'
#SharedGenes Path
sharedGenes_path = matrix_dir + '/sharedGenes.csv'
#Read files into data frame
LGG_Female_df = pd.read_csv(LGG_Female_path, sep='\t')
LGG_Male_df = pd.read_csv(LGG_Male_path, sep='\t')
GBM_Female_df = pd.read_csv(GBM_Female_path, sep='\t')
GBM_Male_df = pd.read_csv(GBM_Male_path, sep='\t')
#Read sharedGenes into data frame
sharedGenes_df = pd.read_csv(sharedGenes_path, index_col = 0)
#Filter by gene ID
LGG_Female_ID = LGG_Female_df['case_id']
LGG_Male_ID = LGG_Male_df['case_id']
GBM_Female_ID = GBM_Female_df['case_id']
GBM_Male_ID =GBM_Male_df['case_id']
#Format LGG_Female ID
LGG_Female_ID_formatted = LGG_Female_ID
for i in range(0,len(LGG_Female_ID)):
   mystring_LGG_Female = LGG_Female_ID[i]
   mystringdots_LGG_Female = mystring_LGG_Female.replace("-", ".")
   print(mystringdots_LGG_Female)
   if mystringdots_LGG_Female[0].isdigit():
      mystringdots_LGG_Female = "X" + mystringdots_LGG_Female
   print(mystringdots_LGG_Female)
   LGG_Female_ID_formatted[i] = mystringdots_LGG_Female
   #Format LGG_Male ID
LGG_Male_ID_formatted = LGG_Male_ID
for i in range(0,len(LGG_Male_ID)):
   mystring_LGG_Male = LGG_Male_ID[i]
   mystringdots_LGG_Male = mystring_LGG_Male.replace("-", ".")
   print(mystringdots_LGG_Male)
   if mystringdots_LGG_Male[0].isdigit():
      mystringdots_LGG_Male = "X" + mystringdots_LGG_Male
   print(mystringdots_LGG_Male)
   LGG_Male_ID_formatted[i] = mystringdots_LGG_Male
   #Format GBM_Female ID
GBM_Female_ID_formatted = GBM_Female_ID
for i in range(0,len(GBM_Female_ID)):
   mystring_GBM_Female = GBM_Female_ID[i]
   mystringdots_GBM_Female = mystring_GBM_Female.replace("-", ".")
   print(mystringdots_GBM_Female)
   if mystringdots_GBM_Female[0].isdigit():
      mystringdots_GBM_Female = "X" + mystringdots_GBM_Female
   print(mystringdots_GBM_Female)
   GBM_Female_ID_formatted[i] = mystringdots_GBM_Female
#Format GBM_Male ID
GBM_Male_ID_formatted = GBM_Male_ID
for i in range(0,len(GBM_Male_ID)):
   mystring_GBM_Male = GBM_Male_ID[i]
   mystringdots_GBM_Male = mystring_GBM_Male.replace("-", ".")
   print(mystringdots_GBM_Male)
   if mystringdots_GBM_Male[0].isdigit():
      mystringdots_GBM_Male = "X" + mystringdots_GBM_Male
   print(mystringdots_GBM_Male)
   GBM_Male_ID_formatted[i] = mystringdots_GBM_Male
#Filter with sharedGenes
sharedGenes_LGG_Female= sharedGenes_df[LGG_Female_ID_formatted.to_list()]
sharedGenes_LGG_Male= sharedGenes_df[LGG_Male_ID_formatted.to_list()]
sharedGenes_GBM_Female= sharedGenes_df[GBM_Female_ID_formatted.to_list()]
sharedGenes_GBM_Male= sharedGenes_df[GBM_Male_ID_formatted.to_list()]
#SAve new files, PATH TO SAVE[]
LGG_Female = matrix_dir + '/sharedGenes_LGG_Female.csv'
LGG_Male = matrix_dir + '/sharedGenes_LGG_Male.csv'
GBM_Female = matrix_dir + '/sharedGenes_GBM_Female.csv'
GBM_Male = matrix_dir + '/sharedGenes_GBM_Male.csv'
#SAVE
sharedGenes_LGG_Female.to_csv(LGG_Female)
sharedGenes_LGG_Male.to_csv(LGG_Male)
sharedGenes_GBM_Female.to_csv(GBM_Female)
sharedGenes_GBM_Male.to_csv(GBM_Male)
print('files have been saved!')