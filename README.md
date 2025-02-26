# Replicating Adebari et al Paper Results

1.  Download the PPI file from GRAND: https://grand.networkmedicine.org/tissues/.
2.  Register for a WebMeV account and download the GBM and LGG data from WebMeV.
3.  Download the curated gene set GMT file from Human MSigDB: https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp.
4.  Download the protein coding genes from GENCODE (basic gene annotation): https://www.gencodegenes.org/human/.
5.  Run **concatenate.py** to concatenate together the GBM and LGG files.
6.  Run **NORMALIZED_DATA_CODE.r** to normalize the data, changing **concatenatedFile** to the path to the concatenated file and **normalizedFile** to the path where you wish to save the normalized file.
7.  Run **Code_to_filter_by_PCG.r** to filter by protein coding genes **mappingFile** to the path to the protein coding genes, **matrixFile** to the path to the concatenated data, and **filteredMatrixFile** to the path where you wish to save the filtered data.
8.  Split the data into male GBM samples, female GBM samples, male LGG samples, and female LGG samples.
9.  Run **runpanda.py** to run PANDA, changing the **localDir** variable to the local directory where the motif, PPI, and expression data are stored.
10.  Run **PercentileMonster.R** to run MONSTER, changing the **localDir** variable to the local directory where the files are stored.
11.  Run **diffanalysis.R** to do pathway analysis changing the **localDir** variable to the local directory where the files are stored and the **pathwayDBFile** to the pathway file path.
12.  Run **Diffanalysis_Visuals.R** to plot the pathways, changing the **localDir** variable to the local directory where the files are stored.
13.  Run **MONSTER_plots.R** changing the **localDir** variable to the local directory where the files are stored.