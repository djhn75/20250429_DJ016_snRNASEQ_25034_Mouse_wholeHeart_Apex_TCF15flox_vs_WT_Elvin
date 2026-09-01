Collection of code and files for the analysis and visualizations presented in the manuscript, related to the TCF15 ChIP-seq data and integration with ATAC-seq:
An endothelial identity program controlled by TCF15 maintains cardiac function
 - TCF15_InfoWriter.py: Script to generate the tables with information on enhancer-gene interaction and information collapsed on gene level.
 - TCF_Exploration.py: Uses the tables from InfoWriter for visualization and additional analyses.
 - MotifEnrichment.py: Run motif enrichment with STREME in TCF15 ChIP-seq peaks with high signal.
 - TCF15_Motif_sites.py: Runs FIMO to find motif hits for the motif found via motif enrichment.
 - HelperScripts: Collection of supplementary functions that are imported by the other scripts.
 - bed_files: Files used and outputted by the above scripts:
   - 0325ATAC_(Diff)Peaks: Used peaks from the ATAC-seq.
   - HUVEC-TCF15_innerJoinMACS3_Signal.bed: Inner join of the TCF15 ChIP-seq peaks from both tags.
   - TCF15Motif_sites.bed: FIMO motif hits for the the motif found via motif enrichment.
   - ext737_CnR_H3K27ac_(Diff)Peaks.bed: Peaks from the H3K27ac ChIP-seq.
   - 0426MarchATACCnR_(Interaction/Gene)Info.txt.gz: Tables that collect the information on interaction/gene level.



