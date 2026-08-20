Collection of code for the analysis and visualizations presented in the manuscript, related to the TCF15 ChIP-seq data and integration with ATAC-seq:
An endothelial identity program controlled by TCF15 maintains cardiac function
 - TCF15_InfoWriter.py: Script to generate the tables with information on enhancer-gene interaction and information collapsed on gene level.
 - TCF_Exploration.py: Uses the tables from InfoWriter for visualization and additional analyses.
 - MotifEnrichment.py: Run motif enrichment with STREME in TCF15 ChIP-seq peaks with high signal.
 - TCF15_Motif_sites.py: Runs FIMO to find motif hits for the motif found via motif enrichment.
 - HelperScripts: Collection of supplementary functions that are imported by the other scripts.



