# Replicating Adebari et al Paper Results

## TCGA
1.  Download the [PPI file](https://granddb.s3.amazonaws.com/tissues/ppi/tissues_ppi.txt) from GRAND.
2.  Register for a WebMeV account and download the GBM and LGG data from WebMeV.
3.  Download the curated gene set GMT file from Human [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp).
4.  Download the protein coding genes from [GENCODE](https://www.gencodegenes.org/human/) (basic gene annotation).
5.  Run **separate_samples.py** to separate data into male and female and find overlapping genes, changing the following:
    -  **gbm_file_path:** The path to the GBM expression data.
    -  **lgg_file_path:** The path to the LGG expression data.
    - **split_file_dir:** The path to the directory where you wish to store the split files.
6.  Run **concatenate.py** to concatenate together the GBM and LGG files.
7.  Run **NORMALIZED_DATA_CODE.r** to normalize the data, changing the following:
    -  **concatenatedFile:** The path to the concatenated file.
    -  **normalizedFile:** The path where you wish to save the normalized file.
8.  Run **Code_to_filter_by_PCG.r** to filter by protein coding gene, changing the following:
    -  **mappingFile:** The path to the protein coding genes.
    -  **matrixFile:** The path to the concatenated data.
    -  **filteredMatrixFile:** The path where you wish to save the filtered data.
9.  Run **updatecode.py** to split the normalized and filtered data into LGG_F, LGG_M, GBM_F, and GBM_M.
    -  **split_file_dir:** The path to the directory where you wish to store the split files.
    -  **matrix_dir:** The path to the file containing the matrix of gene expression levels.
10.  Run **runpanda.py** to run PANDA, changing the following:
     -  **localDir:** The local directory where the motif, PPI, and expression data are stored.
11.  Run **PercentileMonster.R** to run MONSTER, changing the following:
     -  **localDir:** The local directory where the files are stored.
12.  Run **diffanalysis.R** to do pathway analysis changing the following:
     -  **localDir:** The local directory where the files are stored.
     -  **pathwayDBFile:** The pathway file path.
13.  Run **PathwayDifferentialAnalysis.R** to plot the pathways, changing the following:
     -  **sourceDir:** The local directory where the files are stored.
14.  Run **RunBLOBFISHAllTFS.R** to run BLOBFISH and plot the results, changing the following:
     -  **sourceDir:** The local directory where the files are stored.

## REMBRANDT (GSE108474)
1. Download the [clinical](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE108474&format=file&file=GSE108474%5FREMBRANDT%5Fclinical%2Edata%2Etxt%2Egz) and [gene expression](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE108474&format=file&file=GSE108474%5FREMBRANDT%5FGeneExpression%2Etxt%2Egz) data from GEO and unzip them.
2.   Run **preprocessREMBRANDT.r** to normalize the data, separate samples into male and female GBM and LGG, and limit to protein coding genes, changing the following:
    -  **sourceDir:** The path to the REMBRANDT gene expression and clinical data.
    -  **pccFile:** The path to the file containing gene annotations, to obtain protein coding genes
3. Run **runPANDA_REMBRANDT.R** to run PANDA on each subgroup, modifying the following:
    -  **fileDir:** The path to the files REMBRANDT expression files from step 24. Run BLOBFISH on REMBRANDT data.
5. Run **plotBLOBFISH_REMBRANDT_noAR.R** to plot the BLOBFISH results, changing the following:
    -  **pathwayDir:** The path to the TCGA pathway results.
    -  **rembrandtResultDir:** The path to the REMBRANDT BLOBFISH results.
    -  **rembrandtDir:** The path to the REMBRANDT gene expression and clinical data.  
    -  **gmtFilePath:** The path to the file containing mappings from genes to pathways (GMT format).
6. Run **compareBLOBFISH_HIF1A.R** to compare the targeting of specific genes within pathways downstream of HIF1A between male and female GBM in TCGA and REMBRANDT, changing the following:
    -  **pathwayDir:** The path to the TCGA pathway results.
    -  **tcgaDir:** The path to the TCGA gene expression and clinical data.
    -  **rembrandtDir:** The path to the REMBRANDT gene expression and clinical data.
    -  **gmtFilePath:** The path to the file containing mappings from genes to pathways (GMT format).
