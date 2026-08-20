import pandas as pd
from itertools import chain
import json
from pybedtools import BedTool
import HelperScripts.Heatmaps as Heatmaps
import HelperScripts.GOEnrichment as GOEnrichment
import HelperScripts.BasicPlotter as BasicPlotter
import HelperScripts.Bed_Analysis as Bed_Analysis
import HelperScripts.Various as Various
import HelperScripts.GTF_Processing as GTF_Processing


"""Based on the TCF15 ChIP-seq and RNA data for OE and KD and the InfoFiles written from them, explore
the relation of TCF15 binding, ATAC changes and RNA."""

tag = '0426MarchATACCnR'
json_file = "InfoFiles/"+tag+"_Meta_JSON.txt"
info_dict = json.load(open(json_file, 'r'))
plot_path = 'Plots/' + tag
go_dfs_out = "GO_Dfs/" + tag + '_'

gene_df = pd.read_table(info_dict['gene_info'], header=0, sep="\t")
interaction_df = pd.read_table(info_dict['interaction_info'], header=0, sep="\t")
gene_name_map = {x[0]: x[1] for x in gene_df[['Ensembl ID', 'Gene name']].values}

# And the DEGs.
degs = {}
rna_files = Various.fn_patternmatch(info_dict['rna_pattern'])
degs_dir = {}
gene_deg_label = {}  # To have a direct map from gene to which direction of DEG.
full_rna_dfs = {}
for rna_tag in rna_files.values():
    gene_df['RNA_'+rna_tag+' sig'] = [val if pd.isna(val) else '*'*val for val in gene_df['RNA_'+rna_tag+'_DEG']]
    degs['RNA_'+rna_tag] = set(gene_df[gene_df['RNA_'+rna_tag+'_DEG'] == True]['Ensembl ID'].values)
    degs_dir['RNA_'+rna_tag+' down'] = set(gene_df[(gene_df['RNA_'+rna_tag+'_DEG'] == True) & (gene_df['RNA_'+rna_tag+'_log2FC'] < 0)]['Ensembl ID'].values)
    degs_dir['RNA_'+rna_tag + ' up'] = set(gene_df[(gene_df['RNA_'+rna_tag+'_DEG'] == True) & (gene_df['RNA_'+rna_tag+'_log2FC'] > 0)]['Ensembl ID'].values)
    gene_deg_label[rna_tag] = {**{g: 'down DEG' for g in degs_dir['RNA_'+rna_tag+' down']}, **{g: 'up DEG' for g in degs_dir['RNA_'+rna_tag+' up']}}
    full_rna_dfs[rna_tag] = gene_df[~gene_df['RNA_'+rna_tag+"_DEG"].isna()]
rna = 'KD'

# And a set of deregulated vascular genes, split by up and down.
vascular_genes_file = 'VEGFA VEGFR2 signaling_vasculature development_term_UP_DOWN_genes.xlsx'
vascular_genes_df = pd.read_excel(vascular_genes_file, sheet_name='Sheet 1')
vascular_genes = {}
for group_tag, dict_tag in [['DEGs Down', 'Downregulated vascular genes'], ['DEGs Up', 'Upregulated vascular genes']]:
    matched, misses = GTF_Processing.match_gene_identifiers(gene_identifiers=vascular_genes_df[vascular_genes_df['Group'] == group_tag]['Symbols'].values[0].split(','),
                                                            gtf_file=info_dict['annotation'], scopes="symbol, alias, uniprot", fields="ensembl, symbol")
    vascular_genes[dict_tag] = [x['ensembl'] for x in matched.values()]

# -------------------------------------------------------------------------------------------------------------------
# Overlap of DEGs and genes with TCF15 peak
# -------------------------------------------------------------------------------------------------------------------
this_rna_df = full_rna_dfs[rna]
# Also Venns for simplicity.
for nonrna_tag, nonrna_set in [['TCF15 in promoter', set(this_rna_df[this_rna_df['Promoter_TCF15'] > 0]['Ensembl ID'])],  # Supp.4b
                               ['TCF15 in enhancers (pval≤0.05)', set(this_rna_df[this_rna_df['TCF15_inEnhancers p-value'] <= 0.05]['Ensembl ID'])]]:  # Supp.4a
    for rna_dict_tag, rna_dict in [['RNA', {"RNA "+rna+" down": degs_dir['RNA_'+rna+' down'], 'RNA '+rna+' up': degs_dir['RNA_'+rna+' up']}],
                                   ['VEGF', {k: set(val) for k, val in vascular_genes.items()}]]:  # Fig.2i+j
        inter_sets = {**rna_dict, **{nonrna_tag: nonrna_set}}
        BasicPlotter.basic_venn(inter_sets, scaled=False, plot_path=plot_path+rna+ "_" + "DEGs_"+rna_dict_tag+"_"+nonrna_tag, number_size=9, xsize=7, ysize=7, formats=['pdf'], normalize_to=0.15)

# -------------------------------------------------------------------------------------------------------------------
# Location and size of peaks
# -------------------------------------------------------------------------------------------------------------------
tcf15_file = [x for x in info_dict['enhancer_overlap'] if x['tag'] == 'TCF15'][0]['file']
diffatac_file = [x for x in info_dict['enhancer_overlap'] if x['tag'] == 'diffATAC'][0]['file']
# Supp.3f
tcf15_loc = Bed_Analysis.gene_location_bpwise(gtf_file=info_dict['annotation'], bed_dict={"TCF15 peaks": tcf15_file},#, 'ATAC peaks': info_dict['enhancer_file']},
                                              plot_path=plot_path, tss_type='all', palette='tab20', formats=['pdf'])

# Fraction of TCF15 peaks with HUVEC overlap.
# Supp.3d
enhancer_bed = BedTool(info_dict['enhancer_file'])
tcf15_bed = BedTool(tcf15_file)
tcf15_in_huvec = tcf15_bed.intersect(enhancer_bed, u=True)
BasicPlotter.basic_pie({"overlap": len(tcf15_in_huvec), "no overlap": len(tcf15_bed) - len(tcf15_in_huvec)},
                       title='TCF15 peaks overlapping ATAC peaks', palette=['#E1BE6A', '#40B0A6'], numerate=True, legend_perc=True,
                       output_path=plot_path+"ATACov".replace(' ', "_"), legend_title='', formats=['pdf'])

# -------------------------------------------------------------------------------------------------------------------
# Cumulative distribution of log2FC
# -------------------------------------------------------------------------------------------------------------------
rna_run = 'RNA_KD'
this_rna_df = full_rna_dfs[rna_run.split('_')[1]]
plot_dict = {}  # Store all comparisons we want to make.
# Group by significance of TCF15 binding.
# Supp.4e
tcf_sigs = [['TCF15 in promoters', val] for val in this_rna_df[(this_rna_df['Promoter_TCF15'] > 0)][rna_run+'_log2FC'].values]
tcf_sigs += [['TCF15 in enhancers p-val≤0.05', val] for val in this_rna_df[(this_rna_df['TCF15_inEnhancers p-value'] <= 0.05)][rna_run+'_log2FC'].values]
plot_dict["TCF15_sigs"] = {'plot_list': tcf_sigs, 'table_width': 0.55, 'table_x_pos': 1.3, 'xlimit': [-1.5, 1.5], 'add_all': True}

# Always add all genes as comparison.
val_col = rna_run+'_log2FC'
all_genes_list = [['All genes', val] for val in this_rna_df[rna_run+'_log2FC'].values]
for plot_tag, vals in plot_dict.items():
    if vals['add_all']:  # Might be a subset of the full_rna genes, so not relying on the function's add_all.
        main_df = pd.DataFrame(all_genes_list + vals['plot_list'], columns=['Gene set', val_col])
        hue_order = ['All genes'] + sorted(set(main_df['Gene set']) - {"All genes"}) if 'hue_order' not in vals else vals['hue_order']
    else:
        main_df = pd.DataFrame(vals['plot_list'], columns=['Gene set', val_col])
        hue_order = list(main_df['Gene set'].drop_duplicates()) if 'hue_order' not in vals else vals['hue_order']
    BasicPlotter.cumulative_plot(main_df, val_col, 'Gene set', output_path=plot_path + rna_run+'_log2FC_' + plot_tag,
                                     title=None if 'title' not in vals else vals['title'], hue_order=hue_order,
                                     table_width=vals['table_width'], table_x_pos=vals['table_x_pos'],
                                     numerate=True, vertical_line=0, xlimit=vals['xlimit'], formats=['pdf'])

# -------------------------------------------------------------------------------------------------------------------
# GOs
# -------------------------------------------------------------------------------------------------------------------
# And once only for the ones expressed in either condition, meaning that we have in the RNA df.
# Supp.4d
tcf15_rna_sigs = {'TCF15 p-value≤0.05': set(full_rna_dfs['KD'][full_rna_dfs['KD']['TCF15_inEnhancers p-value'] <= 0.05]['Ensembl ID']),}
_ = GOEnrichment.go_enrichment(tcf15_rna_sigs, title_tag="Genes with enrichment of TCF15 in enhancers with KD RNA data",
                               out_tag=plot_path + "EnrichTCF15RNAKD", max_terms='all', font_s=16, numerate=True, wanted_sources=['GO:BP'],
                               background={k: set(full_rna_dfs['KD'][full_rna_dfs['KD']['#Enhancer'] > 0]['Ensembl ID']) for k in tcf15_rna_sigs.keys()},
                               organism='hsapiens', cmap='plasma')

# -------------------------------------------------------------------------------------------------------------------
# Top genes with TCF15 enrichment and predefined sets.
# -------------------------------------------------------------------------------------------------------------------
gene_set_file = 'vascular_and_tgf_markers.xlsx'
gene_set_df = pd.read_excel(gene_set_file, sheet_name=None)
marker_sets = {}
for sheet in gene_set_df:
    matched, misses = GTF_Processing.match_gene_identifiers(gene_identifiers=gene_set_df[sheet]['Unnamed: 0'].values, gtf_file=info_dict['annotation'],
                                                            species='human', scopes="symbol, alias, uniprot", fields="ensembl, symbol")
    marker_sets[sheet + "_All"] = [val['ensembl'] for val in matched.values()]
    matched_spec, misses_spec = GTF_Processing.match_gene_identifiers(gene_identifiers=[x for x in gene_set_df[sheet]['Specific genes'].values if not pd.isna(x)],
                                                                      gtf_file=info_dict['annotation'], species='human', scopes="symbol, alias, uniprot", fields="ensembl, symbol")
    marker_sets[sheet + '_Specific'] = [val['ensembl'] for val in matched_spec.values()]

ranks = 20
s_tag = 'TCF15Peaks'
rna_tag = 'KD'
this_rna_df = full_rna_dfs[rna_tag]
for g_tag, gene_list in [['VascularManual', ['PECAM1', 'CD34', 'KDR', 'CXCR4', 'ICAM1', 'APLN', 'AQP1', 'SOX18']],  # Fig.2k upper
                         ['MesenchymalManual', ['SERPINE1', 'IL6', 'MMP2', 'TNFRSF12A', 'ACTA2', 'TGFB1', 'ITGA5']],  # Fig.2k lower
                         ['SubsetsManual', ['CD34', 'KDR', 'EFCAB14', 'TGFB2', 'RSPO3', 'KDM4C']]] + \
                          list(vascular_genes.items()):  # Supp.5b and Supp.4c
    # Sort again for the cases we use external gene sets.
    if 'Manual' in g_tag:
        to_plot_df = this_rna_df[this_rna_df['Gene name'].isin(gene_list)].set_index('Ensembl ID').sort_values(by=['TCF15_inEnhancers', 'TCF15_inEnhancers p-value'],
                                                                                                               ascending=[False, True])
    else:
        to_plot_df = this_rna_df.set_index('Ensembl ID').loc[[g for g in gene_list if g in full_rna_dfs[rna_tag]['Ensembl ID'].values]].sort_values(by='TCF15_inEnhancers p-value',
                                                                                                   ascending=True)
    to_plot_df = to_plot_df.iloc[:ranks]
    to_plot_df['p-value sig'] = ['*' if val <= 0.05 else '' for val in to_plot_df['TCF15_inEnhancers p-value'].values]
    to_plot_df.columns = [c.replace('_', ' ') for c in to_plot_df]
    to_plot_df['Gene name'] = ['$'+i+'$' for i in to_plot_df['Gene name'].values]

    cmap_cols = {1: {'cols': ['TCF15 inEnhancers p-value'],
                     'cmap': 'viridis',
                     'cbar_label': "p-value"},
                 2: {'cols': ['#Enhancer', 'TCF15 inEnhancers', 'TCF15Motif inEnhancers'],
                     'cmap': 'viridis',
                     'vmin': 0,
                     'cbar_label': "counts"},
                 3: {'cols': ['Promoter TCF15', 'Promoter TCF15Motif'],
                     'vmin': 0,
                     'cmap': 'Greys',
                     'cbar_label': "counts"},
                4: {'cols': ['diffATAC inEnhancers:avg log2FC'],
                    'centre': 0,
                    'cmap': 'bwr',
                    'cbar_label': 'avg log2FC of differential peaks'},
                 7: {'cols': ['RNA ' + rna_tag + ' log2FC'],
                    'centre': 0,
                    'cmap': 'bwr',
                    'cbar_label': "RNA log2FC"},
    }
    Heatmaps.heatmap_cols(to_plot_df, cmap_cols, plot_out=plot_path + g_tag + "_RNA"+rna_tag+'sorted'+s_tag+"Top20",
                          annot_cols=dict({"#Enhancer": "#Enhancer", 'TCF15 inEnhancers': 'TCF15 inEnhancers',
                                           'TCF15Motif inEnhancers': 'TCF15Motif inEnhancers'},
                                          **{k: k for k in list(chain(*[val['cols'] for val in cmap_cols.values()])) if 'Promoter' in k},
                                          **{'RNA '+rna_tag+' log2FC': 'RNA '+rna_tag+' sig', 'TCF15 inEnhancers FDR': 'FDR sig', 'TCF15 inEnhancers p-value': 'p-value sig'}),#, 'RNA '+rna_tag+' log2FC': 'RNA '+rna_tag+' sig'}),
                          annot_s=16, row_label_col="Gene name", class_col=None, x_size=28, y_size=max(5, to_plot_df.shape[0]-8), wspace=0.8, x_rotation=90, formats=['pdf'])

# -------------------------------------------------------------------------------------------------------------------
# Log2FC of the differential peaks and the associated genes
# -------------------------------------------------------------------------------------------------------------------
# Subset for the RNA genes.
rna_gene_interactions = interaction_df[interaction_df['Ensembl ID'].isin(set(full_rna_dfs['KD']['Ensembl ID']))]
# For the replicate information we have to go back to the excel files.
diff_excels = {"H3K27ac": 'ext737-exp1-cnr-punctate-2.xlsx',
               "ATAC": "ext647-exp1-atac-1(1).xlsx"}
diff_cols = {"ATAC": ['kd-siRNA-ctrl_{}'.format(i) for i in range(1, 4)] + ['kd-siRNA-TCF15_{}'.format(i) for i in range(1, 4)],
             "H3K27ac": ['kd-siRNA-ctrl-H3K27ac_{}'.format(i) for i in range(1, 4)] + ['kd-siRNA-TCF15-H3K27ac_{}'.format(i) for i in range(1, 4)]}


def peak_reformat(peak_string):
    if "\t" in peak_string:
        return peak_string.split('\t')[0] + ':' + peak_string.split('\t')[1] + '-' + peak_string.split('\t')[2]
    else:
        return peak_string.split(':')[0] + '\t' + peak_string.split('-')[0].split(':')[1] + '\t' + peak_string.split('-')[1]


# Supp.5a+d
for mark in ['ATAC', "H3K27ac"]:
    diff_interactions = rna_gene_interactions[rna_gene_interactions['ov_diff'+mark] != 0]
    diffpeak_genes = pd.DataFrame(diff_interactions[['peak_id', 'Gene Name']].groupby('peak_id')['Gene Name'].unique())
    diff_bed = BedTool('\n'.join(['\t'.join(x) for x in diff_interactions[['Peak_Chr', 'Peak_Start', "Peak_End"]].astype(str).values]), from_string=True)
    # Get the activities from the Excel file.
    diff_df = pd.read_excel(diff_excels[mark], sheet_name='bigmatrix', header=11, engine='openpyxl')
    diff_df['peak id'] = diff_df['Peak chromosome'] + ':' + diff_df['Peak start'].astype(str) + '-' + diff_df['Peak stop'].astype(str)

    # For H3K27ac we replace the peak name with the overlapping H3K27ac peak.
    if mark == 'H3K27ac':
        excel_bed = BedTool('\n'.join(['\t'.join(x) for x in diff_df[['Peak chromosome', 'Peak start', "Peak stop"]].astype(str).values]), from_string=True)
        overlap_peaks = Bed_Analysis.peaks_peaks_overlap(diff_bed, excel_bed)
        diffpeak_genes.index = [peak_reformat(next(iter(overlap_peaks[peak_reformat(i)]))) for i in diffpeak_genes.index]
        # We have cases where one H3K27ac peak overlapped multiple ATAC peaks.
        diffpeak_genes = pd.DataFrame(diffpeak_genes.groupby(by=diffpeak_genes.index)['Gene Name'].apply(lambda x: list(pd.unique([item for sublist in x for item in sublist]))), columns=['Gene Name'])

    diffpeak_genes['Peak string'] = [i + ': $' + '$, $'.join(val['Gene Name'])+"$" for i, val in diffpeak_genes.iterrows()]
    diff_df = diff_df[diff_df['peak id'].isin(diffpeak_genes.index)].set_index('peak id')
    diffpeak_genes = diffpeak_genes.join(diff_df[diff_cols[mark]])
    diffpeak_genes = diffpeak_genes[['Peak string'] + diff_cols[mark]].set_index('Peak string')
    Heatmaps.clustermap(diffpeak_genes, columns=list(diffpeak_genes.columns), row_column='index', cbar_label="z-score {} signal".format(mark),
                        title="DiffPeaks {}".format(mark), plot_out=plot_path + mark + "_DiffPeaksGenes",
                        cmap='bwr', x_size=18, y_size=12, y_dendro=False, x_dendro=False, main_space=0.5,
                        column_labels=True, row_cluster=True, col_cluster=False, centre=0, tick_size=12,
                        z_score=1, x_rotation=90)


