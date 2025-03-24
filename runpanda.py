import os
from netZooPy.panda.panda import Panda
import pandas as pd
import matplotlib.pyplot as plt
import sys
#Expression data
localDir = null
expression_data_LGG_Female= localDir + '/sharedGenes_LGG_Female.csv'
expression_data_LGG_Male= localDir + '/sharedGenes_LGG_Male.csv'
expression_data_GBM_Female= localDir + '/sharedGenes_GBM_Female.csv'
expression_data_GBM_Male= localDir + '/sharedGenes_GBM_Male.csv'
#Motif data
Female_motif_data     =localDir + '/MotifPriorGencode_p5_female_PARonX.txt'
Male_motif_data = 'localDir + '/MotifPriorGenocode_p5.txt'
#PPI data
ppi_data       =localDir + '/ppi_Gencode_p5.txt'
#Output path
LGG_Female_panda_output   =localDir + '/LGG_Female_Panda_output.txt'
LGG_Male_panda_output   =localDir + '/LGG_Male_Panda_output.txt'
GBM_Female_panda_output   =localDir + '/GBM_Female_Panda_output.txt'
GBM_Male_panda_output   =localDir + '/GBM_Male_Panda_output.txt'
#read expression data
expression_LGG_Female=pd.read_csv(expression_data_LGG_Female,index_col=0)
expression_LGG_Male=pd.read_csv(expression_data_LGG_Male,index_col=0)
expression_GBM_Female=pd.read_csv(expression_data_GBM_Female,index_col=0)
expression_GBM_Male=pd.read_csv(expression_data_GBM_Male,index_col=0)
#read female motif data
Female_motif_data=pd.read_csv(Female_motif_data,sep="\t",header=None)
print(Female_motif_data[0])
#Read male motif data
Male_motif_data=pd.read_csv(Male_motif_data,sep="\t",header=None)
print(Male_motif_data)
#read ppi data
ppi_data=pd.read_csv(ppi_data,sep="\t",header=None)
ppi_data
Female_Motif_mod= Female_motif_data[Female_motif_data[1].isin(expression_LGG_Female.index)]
#call PANDA for LGG_Female
panda_obj_LGG_Female = Panda(expression_LGG_Female, Female_Motif_mod, ppi_data,precision='single', save_tmp=True,save_memory = False, remove_missing=False, keep_expression_matrix = False)
panda_obj_LGG_Female.save_panda_results(LGG_Female_panda_output)
sys.getsizeof(panda_obj_LGG_Female)
#call PANDA for LGG_Male
Male_LGG_Motif_mod= Male_motif_data[Male_motif_data[1].isin(expression_LGG_Male.index)]
panda_obj_LGG_Male = Panda(expression_LGG_Male, Male_LGG_Motif_mod, ppi_data,precision= 'single',  save_tmp=True,save_memory = False, remove_missing=False, keep_expression_matrix = False)
panda_obj_LGG_Male.save_panda_results(LGG_Male_panda_output)
sys.getsizeof(panda_obj_LGG_Male)
#call PANDA for GBM_Female
Female_GBM_Motif_mod = Female_motif_data[Female_motif_data[1].isin(expression_LGG_Female.index)]
panda_obj_GBM_Female = Panda(expression_GBM_Female, Female_GBM_Motif_mod, ppi_data,precision='single', save_tmp=True,save_memory = False, remove_missing=False, keep_expression_matrix = False)
panda_obj_GBM_Female.save_panda_results(GBM_Female_panda_output)
sys.getsizeof(panda_obj_GBM_Female)
#call PANDA for GBM_Male
Male_GBM_Motif_mod= Male_motif_data[Male_motif_data[1].isin(expression_LGG_Male.index)]
panda_obj_GBM_Male = Panda(expression_GBM_Male, Male_GBM_Motif_mod, ppi_data,precision='single', save_tmp=True,save_memory = False, remove_missing=False, keep_expression_matrix = False)
panda_obj_GBM_Male.save_panda_results(GBM_Male_panda_output)
sys.getsizeof(panda_obj_GBM_Male)

