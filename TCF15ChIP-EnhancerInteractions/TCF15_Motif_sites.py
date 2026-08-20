import subprocess
import os
import pandas as pd
from pybedtools import BedTool

"""Get the positions of bHLH motifs within the TCF15 ChIP-seq peaks."""

tcf15_peak_file = "HUVEC-TCF15_innerJoinMACS3.bed"
meme_file = 'TCF15_bHLH_STREME-1.txt'
fimo_out = 'Fimo_runs/TCF15Motif/'
fasta_file = 'hg38.fa'
fimo_src = "meme-5.1.1/src/fimo"
fimo_out_file = fimo_out + "Out_Fimo_TFBSMatrix.txt.gz"  # The name is fixed by the function.

if not os.path.isfile(fimo_out_file):
    # Additionally, call TFBS in the enhancer regions with Fimo, to add those as columns.
    fimo_region_cmd = "python3 FIMO_TFBS_inRegions.py --bed_file " + \
                      tcf15_peak_file + " --PWMs " + meme_file + " --fasta " + fasta_file + \
                      " --fimo_src " + fimo_src + " --out_dir " + fimo_out
    subprocess.call(fimo_region_cmd, shell=True)

# We don't want the counts per TCF15 peak, but the location of the predicted TFBS themselves.
tfbs_df = pd.read_table(fimo_out + 'Out_Fimo.tsv.gz', sep='\t', header=0)
tfbs_df['chr'] = [x.split('::')[1].split(':')[0] for x in tfbs_df['sequence_name'].values]
tfbs_df['tfbs_start'] = [int(val['sequence_name'].split(':')[3].split('-')[0]) + val['start'] for val in tfbs_df.to_dict(orient='records')]
tfbs_df['tfbs_end'] = [int(val['sequence_name'].split(':')[3].split('-')[0]) + val['stop'] for val in tfbs_df.to_dict(orient='records')]
tfbs_sites = BedTool('\n'.join(['\t'.join([str(y) for y in x]) for x in tfbs_df[['chr', 'tfbs_start', 'tfbs_end']].values]), from_string=True).sort().merge()
open(fimo_out + '/TCF15Motif_sites.bed', 'w').write(str(tfbs_sites))



