from pybedtools import BedTool
import numpy as np
import HelperScripts.Various as Various
import HelperScripts.Bed_Analysis as Bed_Analysis
import HelperScripts.MEME_Shorts as MEME_Shorts

"""Collect the motif enrichment runs here. We also need preprocessing of the files to get the summit for 
the TCF15 ChIP-seq."""
sequence_file = 'hg38.fa'
streme_out = 'STREME_runs/'
sea_out = 'SEA_runs/'
meme_file = 'PWMs/2.2/Jaspar_Hocomoco_Kellis_human_meme.txt'

# ------------------------------------------------------------------------------------------------------------
# STREME in TCF15 ChIP-seq peaks
# ------------------------------------------------------------------------------------------------------------
# First get the summit position of the three narrowpeak files along with their signalValue.
narrowpeak_files = Various.fn_patternmatch("HUVEC-TCF15-*_fdr0.1_peaks.narrowPeak")
narrowpeak_coll = {}
for file, tag in narrowpeak_files.items():
    narrow = BedTool(file)
    narrowpeak_coll[file] = {'\t'.join(x.fields[:3]): {"summit": int(x.fields[-1]), 'signal': float(x.fields[6])} for x in narrow}

p_tag = 'InnerJoin'
p_file = 'HUVEC-TCF15_innerJoinMACS3.bed'
p_bed = BedTool(p_file)
p_summits = {'\t'.join(x.fields[:3]): {'summit': None, 'signal': 0, 'all_signals': []} for x in p_bed}
for narrow in narrowpeak_files:
    narrow_map = Bed_Analysis.peaks_peaks_overlap(p_bed, narrow)
    for p_hit, n_hits in narrow_map.items():
        if n_hits:
            for n_hit in n_hits:  # Take the summit with the highest signal.
                if p_summits[p_hit]['signal'] < narrowpeak_coll[narrow][n_hit]['signal']:
                    p_summits[p_hit]['signal'] = narrowpeak_coll[narrow][n_hit]['signal']
                    p_summits[p_hit]['summit'] = int(n_hit.split('\t')[1]) + narrowpeak_coll[narrow][n_hit]['summit']
                p_summits[p_hit]['all_signals'].append(narrowpeak_coll[narrow][n_hit]['signal'])

# After getting the summit with the highest signal, take the top peaks, based on the average signal of the
# overlapping peaks.
top_peaks = sorted([x for x in p_summits], key=lambda x: np.mean(p_summits[x]['all_signals']), reverse=True)[:1000]
top_peaks_summits = []
outsiders = []
for peak in top_peaks:
    summit_loc = p_summits[peak]['summit']
    if summit_loc > int(peak.split('\t')[2]) or summit_loc < int(peak.split('\t')[1]):
        outsiders.append([peak, summit_loc])
    top_peaks_summits.append([peak.split('\t')[0], str(summit_loc - 100), str(summit_loc + 100)])
print(p_tag, 'outliers', len(outsiders))

summits_bed = BedTool('\n'.join(['\t'.join(x) for x in top_peaks_summits]), from_string=True)

MEME_Shorts.streme(foreground=summits_bed, sequence=sequence_file,
                   file_tag='STREME_MACS3'+p_tag, out_dir=streme_out)

