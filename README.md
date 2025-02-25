# Replicating Adebari et al Paper Results

1.  Download the PPI file from GRAND: https://grand.networkmedicine.org/tissues/.
2.  Register for a WebMeV account and download the GBM and LGG data from WebMeV.
3. Download the curated gene set GMT file from Human MSigDB: https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp.
4.  Run concatenate.py to concatenate together the GBM and LGG files.
4.  Run runpanda.py to run PANDA, changing the **localDir** variable to the local directory where the motif, PPI, and expression data are stored.
4.  Run PercentileMonster.R to run MONSTER, changing the **localDir** variable to the local directory where the files are stored.
4.  Run diffanalysis.R to do pathway analysis changing the **localDir** variable to the local directory where the files are stored and the **pathwayDBFile** to the pathway file path.
4.  Run Diffanalysis_Visuals.R to plot the pathways, changing the **localDir** variable to the local directory where the files are stored.
5.  Run MONSTER_plots.R changing the **localDir** variable to the local directory where the files are stored.