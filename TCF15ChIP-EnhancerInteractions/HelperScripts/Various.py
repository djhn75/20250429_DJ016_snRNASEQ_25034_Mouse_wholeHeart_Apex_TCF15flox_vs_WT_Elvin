from collections import Counter
import pandas as pd
from multiprocess import Pool
import scipy.stats
import statsmodels.stats.multitest
import os
import fnmatch
import re


"""Collection of functions that didn't fit well anywhere else. Also includes functions which are
required by multiple other modules."""


def sanitize_path(path_string):
    """Function to remove unwanted characters from a file path before saving. Put here since we need it for a lot
    of functions and eases adjustment e.g. for OS."""
    return re.sub(r'[^A-Za-z0-9._\-\\\\////]+', '', path_string)
    # return path_string.replace(' ', '').replace(':', '').replace("°", "Degree").replace()


def fn_patternmatch(pattern):
    """
    Grabs all files in the file system that match the pattern and returns a dictionary with {file: wildcard}.
    Only works if the wildcard is in the file name and not in a directory name.
    Files that start with ._ (e.g. temporary or system files) will be skipped.
    E.g. for a directory BirdCollection/ that contains Bird_Kakapo.txt and Bird_Kea.txt:
    fn_patternmatch("BirdCollection/*.txt" = {"BirdCollection/Bird_Kakapo.txt": "Kakao", "BirdCollection/Bird_Kea.txt": "Kea"}
     """
    parent_folder = '/'.join(pattern.split('/')[:-1]) + '/'
    children_pattern = pattern.split('/')[-1]
    re_pattern = re.compile(children_pattern.replace('*', '(.*)'))
    matched_files = {parent_folder + x: re_pattern.search(x).group(1)
                     for x in os.listdir(parent_folder)
                     if fnmatch.fnmatch(x, children_pattern) and not x.startswith('._')}
    return matched_files


def cre_fisher(args):
    """Separate function to enable parallelization."""
    gene, fish_table = args
    _, pval = scipy.stats.fisher_exact(fish_table, alternative='greater')
    return [gene, pval]


def gene_cre_overlap_fisher(interaction_df, overlap_col, peak_id_col='peak_id', ncores=1):
    """CARE: the input df is quite specific.
    Takes a pandas df of interactions (gene in Ensembl ID and a peak identifier in peak_id_col) and tests for each gene
    in the Df if its CREs more often have a hit in overlap_col (e.g. ChIP-seq peak overlap) than compared to the whole
    set of CREs that form interactions."""
    inter_genes = list(set(interaction_df['Ensembl ID']))
    inter_cres = len(set(interaction_df[peak_id_col]))
    gene_num_cres = Counter(interaction_df['Ensembl ID'])
    hits_only = interaction_df[interaction_df[overlap_col] > 0]
    cres_w_hits = len(set(hits_only[peak_id_col]))
    gene_cre_hits = Counter(hits_only['Ensembl ID'])
    fisher_tables = []
    for gene in inter_genes:
        fisher_tables.append([[gene_cre_hits[gene], gene_num_cres[gene] - gene_cre_hits[gene]],
                              [cres_w_hits - gene_cre_hits[gene], inter_cres - cres_w_hits - (gene_num_cres[gene] - gene_cre_hits[gene])]])
    
    process_pool = Pool(processes=ncores)
    fish_pool = process_pool.map(cre_fisher, [[g, fisher_tables[i]] for i, g in enumerate(inter_genes)])
    process_pool.close()
    fisher_df = pd.DataFrame(fish_pool, columns=['Ensembl ID', overlap_col+' p-value']).set_index('Ensembl ID')
    fisher_df[overlap_col+' FDR'] = statsmodels.stats.multitest.fdrcorrection(fisher_df[overlap_col+' p-value'],
                                                                              alpha=0.05, method='indep', is_sorted=False)[1]
    return fisher_df


