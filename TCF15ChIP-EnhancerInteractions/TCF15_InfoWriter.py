from pybedtools import BedTool
import gzip
import pandas as pd
import os
import json
import subprocess
from collections import Counter
import numpy as np
from timeit import default_timer as clock
import HelperScripts.Bed_Analysis as Bed_Analysis
import HelperScripts.Various as Various

"""Create information files that aggregate the different files we have, one for the genes and on the level 
of enhancer-gene interactions."""

# First we need to give the paths to all the different files.
info_dict = {'out_tag': '0426MarchATACCnR',
             'interaction_folder': "0325ATAC_ABC/",
             'abc_cols': ['Contact', 'TSS-dist', 'ABC-Score'],
             'enhancer_file': "0325ATAC_Peaks.bed",
             'rna_pattern': "TCF15_*_RNA_formatted.txt",
             'out_path': "InfoFiles/",
             'enhancer_in_promoter': True,
             'enhancer_overlap': [{'tag': 'TCF15', 'file': 'HUVEC-TCF15_innerJoinMACS3_Signal.bed', 'same_peak': False, 'fetch': 'Signal'},
                                  {'tag': 'TCF15Motif', 'file': 'TCF15Motif_sites.bed', 'same_peak': False, 'fetch': None},
                                  {'tag': 'diffATAC', 'file': '0325ATAC_DiffPeaks.bed', 'same_peak': True, 'fetch': 'log2FC', 'no_fisher': True},
                                  {'tag': 'H3K27ac', 'file': 'ext737_CnR_H3K27ac_Peaks.bed', 'same_peak': False, 'fetch': None},
                                  {'tag': 'diffH3K27ac', 'file': 'ext737_CnR_H3K27ac_DiffPeaks.bed', 'same_peak': False, 'fetch': 'log2FC', 'no_fisher': True},
                                  ],
             'promoter_overlap': [{'tag': 'TCF15', 'file': 'HUVEC-TCF15_innerJoinMACS3_Signal.bed', 'fetch': 'Signal'},
                                  {'tag': 'TCF15Motif', 'file': 'TCF15Motif_sites.bed', 'fetch': None},
                                  {'tag': 'diffATAC', 'file': '0325ATAC_DiffPeaks.bed', 'fetch': 'log2FC'},
                                  {'tag': 'H3K27ac',
                                   'file': 'ext737_CnR_H3K27ac_Peaks.bed',
                                   'same_peak': False, 'fetch': None},
                                  {'tag': 'diffH3K27ac',
                                   'file': 'ext737_CnR_H3K27ac_DiffPeaks.bed',
                                   'same_peak': False, 'fetch': 'log2FC'},
                                  ],
             'genebody_overlap': [],
             'annotation': "gencode.v38.annotation.gtf.gz",
             }

info_dict['interaction_info'] = info_dict['out_path'] + info_dict['out_tag'] + '_InteractionInfo.txt.gz'
info_dict['gene_info'] = info_dict['out_path'] + info_dict['out_tag'] + '_GeneInfo.txt.gz'

# Write the dictionary to a file, to later be able to reconstruct which files came from where, and to ease the
# follow-up plotting functions.
open(info_dict['out_path'] + info_dict['out_tag'] + '_Meta_JSON.txt', 'w').write(json.dumps(info_dict))

gene_name_map = {g.strip().split('\t')[8].split('gene_id "')[-1].split('";')[0].split('.')[0]:
                 g.strip().split('\t')[8].split('gene_name "')[-1].split('";')[0] for g in gzip.open(info_dict['annotation'], 'rt').readlines() if not g.startswith('#') and g.split('\t')[2] == 'gene'}

interaction_files = {info_dict['interaction_folder'] + '/' + x: x.split('ABCpp_scoredInteractions_')[-1].split('.')[0]
                     for x in os.listdir(info_dict['interaction_folder']) if "ABCpp_scoredInteractions" in x and x.endswith(".gz")}

# For each enhancer get the information we already have.
enhancer_bed = BedTool(info_dict['enhancer_file'])
enhancer_dict = {'\t'.join([str(y) for y in k]): {} for k, vals in
                 pd.read_table(info_dict['enhancer_file'], header=None, dtype=str, index_col=[0, 1, 2]).to_dict(orient="index").items() if not k[0].startswith('#')}

# ---------------------------------------------------------------------------------------------------
# Interaction-Info
# ---------------------------------------------------------------------------------------------------
# Define which columns we want to have in the output file, and append it when adding new information.
# Add the peak_id entry.
if 'peak_id' not in info_dict or not info_dict['peak_id']:
    for enh in enhancer_dict:
        enhancer_dict[enh]['peak_id'] = enh.split('\t')[0] + ':' + enh.split('\t')[1] + '-' + enh.split('\t')[2]
else:
    for enh in enhancer_dict:
        enhancer_dict[enh]['peak_id'] = enhancer_dict[enh][info_dict['peak_id']]

# Define promoter regions around the 5' TSS of all annotated genes in the gtf-file.
if 'enhancer_in_promoter' in info_dict and info_dict['enhancer_in_promoter']:
    peak_ovs, _ = Bed_Analysis.peaks_promoter_overlap(peak_file=info_dict['enhancer_file'], gtf_file=info_dict['annotation'], tss_type='all')
    for peak, genes in peak_ovs.items():
        enhancer_dict[peak]['ov_Promoter'] = genes

if 'enhancer_overlap' in info_dict:
    for peak_vals in info_dict['enhancer_overlap']:
        peak_ovs = Bed_Analysis.peaks_peaks_overlap(peak_file=info_dict['enhancer_file'],
                                                    other_peak_file=peak_vals['file'])
        for peak, hits in peak_ovs.items():
            enhancer_dict[peak]['ov_'+peak_vals['tag']] = len(hits)
        if peak_vals['fetch']:
            peak_hits, hit_cols = Bed_Analysis.peaks_fetch_col(base_regions=info_dict['enhancer_file'],
                                                               pattern=peak_vals['file'], same_peaks=peak_vals['same_peak'], fetch_col=peak_vals['fetch'])
            for p_hit, hit in peak_hits.items():
                enhancer_dict[p_hit]['ov_'+peak_vals['tag']+':'+peak_vals['fetch']] = hit[hit_cols[0]]  # We only have one matching file.

# We collect the enhancer-gene interactions from conditions, and store their ABC-scores.
interaction_set = {}
start = clock()
for file in interaction_files:
    with gzip.open(file, 'rt') as inter_open:
        inter_head = {x: i for i, x in enumerate(inter_open.readline().strip().split('\t'))}
        for entry in inter_open:
            entry = entry.strip().split('\t')
            this_enhancer = '\t'.join([entry[inter_head[c]] for c in ['#chr', 'start', 'end']])
            this_inter = entry[inter_head['Ensembl ID']].split('.')[0] + '#' + entry[inter_head['Gene Name']] + "#" + this_enhancer
            if this_inter not in interaction_set:
                interaction_set[this_inter] = {}
            interaction_set[this_inter][interaction_files[file] + ' ABC-Score'] = entry[inter_head['ABC-Score']]
            if 'abc_cols' in info_dict and info_dict['abc_cols']:
                for col in info_dict['abc_cols']:
                    interaction_set[this_inter][col] = entry[inter_head[col]]
print(clock() - start)

# The output writing is admittedly ugly, but storing everything on interaction level right away is highly redundant.
# And then having to look up the enhancers iteratively wouldn't be efficient with datastructures like pandas DFs.
abc_header = ['Peak_Chr', 'Peak_Start', 'Peak_End', 'Ensembl ID', 'Gene Name']
abc_cols = [x for x in info_dict['abc_cols'] if x != 'ABC-Score'] + [c+" ABC-Score" for c in interaction_files.values()]
abc_header += abc_cols
all_enh_cols = list(sorted(set(next(iter(enhancer_dict.values())).keys())))
abc_header += all_enh_cols
with open(info_dict['interaction_info'].replace('.gz', ''), 'w') as output:
    output.write('\t'.join(list(abc_header)) + '\n')
    for inter, inter_vals in interaction_set.items():
        if inter.split('#')[2].startswith('chr'):
            enh = inter.split('#')[2]
        else:
            enh = 'chr' + inter.split('#')[2]
        if enh in enhancer_dict:
            output.write(enh + '\t' + inter.split('#')[0] + '\t' +
                         gene_name_map[inter.split('#')[0].split('.')[0]])
            for col in abc_cols:
                output.write('\t' + ('0' if col not in inter_vals else inter_vals[col]))
            enh_vals = enhancer_dict[enh]
            for col in all_enh_cols:
                output.write('\t' + (str(enh_vals[col]) if type(enh_vals[col]) not in [set, list] else ','.join(sorted(enh_vals[col]))))
            output.write('\n')


# ---------------------------------------------------------------------------------------------------
# Gene-Info
# ---------------------------------------------------------------------------------------------------
# Prepare the datastructure for aggregating interaction information.
# Write the output for all genes in the annotation.
gene_df = pd.DataFrame()
gene_df['Ensembl ID'] = list(gene_name_map.keys())
gene_df['Gene name'] = [gene_name_map[g] for g in gene_df['Ensembl ID']]
# Switch to either open or expressed
# Add the RNA tables as join.
rna_files = Various.fn_patternmatch(info_dict['rna_pattern'])
for rna_file, rna_tag in rna_files.items():
    rna_df = pd.read_table(rna_file, sep='\t', header=0).set_index('Ensembl ID').drop("Gene name", axis=1)
    rna_df.columns = ['RNA_'+rna_tag+'_'+c for c in rna_df.columns]
    gene_df = gene_df.join(rna_df, on='Ensembl ID', how='left')

for peak_vals in info_dict['promoter_overlap']:
    _, gene_hits = Bed_Analysis.peaks_promoter_overlap(peak_file=peak_vals['file'], gtf_file=info_dict['annotation'],
                                                       tss_type='all')
    gene_df['Promoter_'+peak_vals['tag']] = [len(gene_hits[g]) for g in gene_df['Ensembl ID']]
    if peak_vals['fetch']:
        gene_hits, hit_cols = Bed_Analysis.promoter_fetch_col(pattern=peak_vals['file'], gtf_file=info_dict['annotation'],
                                                              fetch_col=peak_vals['fetch'])
        gene_df['Promoter_'+peak_vals['tag']+':'+peak_vals['fetch']] = [gene_hits[g][hit_cols[0]] for g in gene_df['Ensembl ID'].values]

for peak_vals in info_dict['genebody_overlap']:
    gb_fracs = Bed_Analysis.peaks_genebody_overlap(peak_file=peak_vals['file'], gtf_file=info_dict['annotation'])
    gene_df['Fraction gene body with ' + peak_vals['tag']] = [gb_fracs[g] for g in gene_df['Ensembl ID']]

# Collect the interaction information on gene level.
interaction_df = pd.read_table(info_dict['interaction_info'].replace('.gz', ''), header=0, sep='\t')
enhancer_counter = Counter(interaction_df['Ensembl ID'])
gene_df['#Enhancer'] = [0 if g not in enhancer_counter else enhancer_counter[g] for g in gene_df['Ensembl ID']]
for peak_vals in info_dict['enhancer_overlap']:
    summed_counts = interaction_df[['Ensembl ID', 'ov_' + peak_vals['tag']]].groupby('Ensembl ID').sum().to_dict(orient='index')
    gene_df[peak_vals['tag']+"_inEnhancers"] = [0 if g not in summed_counts else summed_counts[g]['ov_' + peak_vals['tag']] for g in gene_df['Ensembl ID']]
    if peak_vals['fetch']:
        avg_values = interaction_df[['Ensembl ID', 'ov_' + peak_vals['tag']+":"+peak_vals['fetch']]].dropna().groupby('Ensembl ID').mean().to_dict(orient='index')
        gene_df[peak_vals['tag']+'_inEnhancers:avg_'+peak_vals['fetch']] = [0 if g not in avg_values else avg_values[g]['ov_' + peak_vals['tag']+":"+peak_vals['fetch']]
                                             for g in gene_df['Ensembl ID'].values]

# Add a Fisher test to the enhancer_overlap data.
for peak_vals in info_dict['enhancer_overlap']:
    if peak_vals['fetch']:  # Then we need to translate the values to ov or not.
        interaction_df['ov_'+peak_vals['tag']] = [~np.isnan(x)*1 for x in interaction_df['ov_'+peak_vals['tag']+':'+peak_vals['fetch']].values]
    if 'no_fisher' not in peak_vals or not peak_vals['no_fisher']:
        fisher_df = Various.gene_cre_overlap_fisher(interaction_df, overlap_col='ov_'+peak_vals['tag']).rename({'ov_'+peak_vals['tag']+" p-value": peak_vals['tag']+'_inEnhancers p-value',
                                                                                                   'ov_'+peak_vals['tag']+" FDR": peak_vals['tag']+'_inEnhancers FDR'}, axis=1)
        gene_df = gene_df.join(fisher_df, how='left', on='Ensembl ID')
        gene_df.loc[gene_df[peak_vals['tag']+"_inEnhancers p-value"].isna(), [peak_vals['tag']+'_inEnhancers p-value', peak_vals['tag']+'_inEnhancers FDR']] = 1

gene_df.to_csv(info_dict['gene_info'].replace('.gz', ''), sep='\t', header=True, index=False)

# Gzip the output files.
subprocess.call("gzip -f " + info_dict['interaction_info'].replace('.gz', ''), shell=True)
subprocess.call("gzip -f " + info_dict['gene_info'].replace('.gz', ''), shell=True)

