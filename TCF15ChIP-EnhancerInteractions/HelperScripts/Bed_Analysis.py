import numpy as np
from pybedtools import BedTool
import pandas as pd
import gzip
import GTF_Processing
import BasicPlotter
import Various

"""Collection of functions related to analysis of bed-files, their location or intersection with other files."""


def gene_location_bpwise(bed_dict, gtf_file, plot_path, tss_type='5', external_bed={}, palette='tab20',
                         formats=['pdf']):
    """
    Based on a gtf file builds bed-objects for Promoter (±200bp), Exons, UTRs and Introns, and then counts how many
    of the bp in the bed file(s) are located within those annotations and which are intergenic. All gene features are
    exclusive, overlaps are removed. Introns are gene bodies subtracted by all other features. The bed_dict can also
    be a list of bed files, to omit recreating the gtf-annotations each time. Creates a pie chart with
    the percentages in the given path. Also returns a dictionary with the bp-location of each bed-region and the
    total overlaps.

    Args:
        bed_dict: A dictionary with {title tag: bed-file or BedTools object}
        gtf_file: gtf-file in GENCODE's format, can be gzipped.
        tss_type: What to consider as promoter of genes, either '5' to use only the most 5' TSS, or 'all' to consider
            ±200bp around all annotated TSS of a gene.
        external_bed: An additional dictionary of bed-file or BedTools object which will be added as category for the
            intersection. This will be considered as highest priority, meaning the regions in there are removed from the
            gene-related features, and a bp overlapping external_bed will not be counted anywhere else. Multiple external
            bed-regions shouldn't overlap, that causes undefined outcomes.

    Returns:
        tuple:
            - **regions_locs**: For each entry in the bed_dict a dictionary with the bp-wise locations for each individual region.
            - **total_locs**: For each entry in the bed_dict the overall overlap of base pairs with each genomic feature.
    """
    if gtf_file.endswith('.gz'):
        bed_annotations = {'Promoter': GTF_Processing.gene_window_bed(gtf_file, 200, tss_type=tss_type).sort().merge(),
                           'Exons': BedTool('\n'.join(['\t'.join([x.strip().split('\t')[c] for c in [0, 3, 4]]) for x in
                                                       gzip.open(gtf_file, 'rt').readlines() if
                                                       not x.startswith('#') and x.split('\t')[2] == 'exon']),
                                            from_string=True).sort().merge(),
                           'UTR': BedTool('\n'.join(['\t'.join([x.strip().split('\t')[c] for c in [0, 3, 4]]) for x in
                                                     gzip.open(gtf_file, 'rt').readlines() if
                                                     not x.startswith('#') and x.split('\t')[2] == 'UTR']),
                                          from_string=True).sort().merge()}
    else:
        bed_annotations = {'Promoter': GTF_Processing.gene_window_bed(gtf_file, 200, tss_type=tss_type).sort().merge(),
                           'Exons': BedTool('\n'.join(['\t'.join([x.strip().split('\t')[c] for c in [0, 3, 4]]) for x in
                                                       open(gtf_file).readlines() if
                                                       not x.startswith('#') and x.split('\t')[2] == 'exon']),
                                            from_string=True).sort().merge(),
                           'UTR': BedTool('\n'.join(['\t'.join([x.strip().split('\t')[c] for c in [0, 3, 4]]) for x in
                                                     open(gtf_file).readlines() if
                                                     not x.startswith('#') and x.split('\t')[2] == 'UTR']),
                                          from_string=True).sort().merge()}
    # Remove the Promoter regions and UTRs from Exons and Promoter from UTRs, want to have exclusive annotations.
    bed_annotations['Exons'] = bed_annotations['Exons'].subtract(bed_annotations['Promoter']).subtract(
        bed_annotations['UTR'])
    bed_annotations['UTR'] = bed_annotations['UTR'].subtract(bed_annotations['Promoter'])

    introns_bed = BedTool(gtf_file).sort().merge()
    for remove in ['Promoter', 'Exons', 'UTR']:
        introns_bed = introns_bed.subtract(bed_annotations[remove])
    bed_annotations['Introns'] = introns_bed

    # Add potential external bed files and subtract it from all other regions.
    if external_bed:
        for external in external_bed:
            bed_annotations[external] = BedTool(external_bed[external]).sort().merge()
            for anno in bed_annotations:
                if anno not in external_bed:
                    bed_annotations[anno] = bed_annotations[anno].subtract(bed_annotations[external])

    regions_locs = {}
    total_locs = {}
    for tag, this_bed in bed_dict.items():
        print(tag)
        if len(BedTool(this_bed)) == 0:
            print("Empty bed", tag)
            continue
        # Merge the bedfile to have unique bp, but keep track of what got merged and assemble it back later.
        bed = BedTool(this_bed).sort().merge(c=[1, 2, 3], o=['collapse'] * 3)
        org_beds = {'\t'.join(x.fields[:3]): [] for x in bed}
        for site in bed:
            org_beds['\t'.join(site.fields[:3])] += [
                '\t'.join([site.fields[3].split(',')[i], site[4].split(',')[i], site[5].split(',')[i]]) for i in
                range(site.fields[4].count(',') + 1)]
        bed_inter = {'\t'.join(x.fields[:3]): {a: 0 for a in bed_annotations.keys()} for x in bed}
        # Very broad marks rarely intersect only one annotation, so we count each bp.
        overall_inter = {a: 0 for a in list(bed_annotations.keys()) + ['Intergenic']}

        for annot, annot_bed in bed_annotations.items():
            intersection = bed.intersect(annot_bed, wo=True)
            for inter in intersection:
                bed_inter['\t'.join(inter.fields[:3])][annot] += int(inter.fields[-1])
                overall_inter[annot] += int(inter.fields[-1])

        # Intergenic is every base that was not assigned to any of the other annotations yet.
        for entry in bed_inter:
            missing_bp = abs(int(entry.split('\t')[2]) - int(entry.split('\t')[1])) - sum(bed_inter[entry].values())
            if missing_bp < 0:
                print("WARNING: There's a negative number of missing bp!")
            overall_inter['Intergenic'] += missing_bp
            bed_inter[entry]['Intergenic'] = missing_bp

        locations = {k: 0 for k in overall_inter.keys()}
        for annot in set(overall_inter.keys()):
            max_inter = len([x for x in bed_inter if max(bed_inter[x], key=bed_inter[x].get) == annot])
            print(annot, round(overall_inter[annot] / sum(overall_inter.values()) * 100, 2), 'max', max_inter,
                  round(max_inter / len(bed) * 100, 2))
            locations[annot] = round(overall_inter[annot] / sum(overall_inter.values()) * 100, 2)

        loc_df = pd.DataFrame([[annot, locations[annot]] for annot in overall_inter.keys()],
                              columns=['Location', 'Overlap']).set_index('Location')
        BasicPlotter.basic_pie(plot_df=loc_df, title=tag + '\n#' + str(len(BedTool(this_bed))), palette=palette,
                               numerate=False, legend_perc=True, formats=formats, legend_title='',
                               output_path=plot_path + (tag + "_GeneFeatureLocation_bpwise").replace(" ", '_'))

        # Now map the potentially merged regions back to all its original regions.
        org_inters = {}
        for inter, locs in bed_inter.items():
            for sub_i in org_beds[inter]:
                org_inters[sub_i] = locs
        regions_locs[tag] = org_inters
        total_locs[tag] = locations

    return regions_locs, total_locs


def peaks_peaks_overlap(peak_file, other_peak_file):
    """
    Intersects two bed-files or BedTools object to return a dictionary with {chr\tstart\tend: {other peaks}}.

    Args:
        peak_file: Path to a bed-file or BedTools object of the regions that will be the keys in the mapping dictionary.
        other_peak_file: Path or BedTools object to be mapped to peak_file and listed as values in the dictionary.

    Returns:
        dict:
            - **peak_dict**: A dictionary with {chr\tstart\tend: {other peaks}}. Values are empty if there was no intersection.
    """
    peak_dict = {'\t'.join(x.fields[:3]): set() for x in BedTool(peak_file)}
    peaks_inter = BedTool(peak_file).intersect(BedTool(other_peak_file), wo=True)
    other_start = len(BedTool(peak_file)[0].fields)
    for inter in peaks_inter:
        peak_dict['\t'.join(inter.fields[:3])].add('\t'.join(inter.fields[other_start:other_start+3]))
    return peak_dict


def peaks_promoter_overlap(peak_file, gtf_file, tss_type='all', gene_set=(), extend=200):
    """
    Based on a bed-file path or BedTools object returns a dictionary with {chr\tstart\tend: {genes whose promoter overlap}}
    and one with {gene: {peaks}}.

    Args:
        peak_file: Path to a bed-file or BedTools object of the regions that will be used for the intersection.
        gtf_file: gtf-file in GENCODE's format, can be gzipped.
        tss_type: "5" to do the overlap only for the 5' TSS or "all" to do it for all unique TSS of all transcripts in the gtf-file.
        gene_set: Set of Ensembl IDs or gene names or mix of both to limit the output to. If empty, return for all
            genes in the annotation.
        extend: Number of base pairs to extend the TSS in each direction. 200 means a window of size 401.

    Returns:
        tuple:
            - **peak_dict**: Dictionary mapping peaks to gene promoters {chr\tstart\tend: {genes whose promoter overlap}}.
            - **gene_dict**: Dictionary doing the mapping the other way around with {gene: {peaks}}.
    """
    promoter_bed = GTF_Processing.gene_window_bed(gtf_file=gtf_file, extend=extend, tss_type=tss_type, merge=True,
                                                  gene_set=gene_set)
    gene_dict = {x.fields[3]: set() for x in promoter_bed}
    peak_file = BedTool(peak_file)
    peak_dict = {'\t'.join(x.fields[:3]): set() for x in peak_file}
    promoter_inter = peak_file.intersect(promoter_bed, wo=True)
    for inter in promoter_inter:
        peak_dict['\t'.join(inter.fields[:3])].add(inter.fields[-4])
        gene_dict[inter.fields[-4]].add('\t'.join(inter.fields[:3]))
    return peak_dict, gene_dict


def peaks_fetch_col(base_regions, pattern, same_peaks=False, fetch_col='log2FC'):
    """
    Take a bed-file or BedTools object with regions as peaks and intersect it with other peak files
    defined with a filesystem pattern to get their fetch_col values in the base_regions. Useful for example, if one has
    multiple differential peak calls and want to map them to a base set of regions. If multiple regions overlap a base
    region, their average signal is taken. Entries of base peaks without overlap will be filled with NaN.

    Args:
        base_regions: BedTool's object or path to a bed file with the regions on which the intersection is centred on.
        pattern: File path pattern e.g., DiffPeaks/DiffBind_*.bed, or the path to just one individual file.
            All files matching that pattern will be used for the intersection. The string at the asterisk will be
            used as identifier. E.g., DiffPeaks/DiffBind_Macrophages.bed will be identified as Macrophages.
            The files need to have a header starting with # with fetch_col as column name. If it's not a path pattern
            but an individual file, the suffix after the last '/' is used as identifier.
        same_peaks: Boolean if the base regions are the same regions as the ones found in the pattern files.
        fetch_col: The column name in the files found by pattern to identify from which column the value should be taken.

    Returns:
        tuple:
            - **fill_dict**: Nested dictionary of the base_regions as keys and the average of the fetch_col in each of the matched files.
            - **matched_files**: List of the identifiers derived from the files matching the pattern.
    """
    fill_bed = BedTool(base_regions)
    fill_dict = {'\t'.join(x.fields[:3]): {} for x in fill_bed}
    if '*' in pattern:
        comparison_files = Various.fn_patternmatch(pattern)
    else:  # Assume we have just one file.
        comparison_files = {pattern: pattern.split('/')[-1]}
    print(comparison_files)
    # Store whether an enhancer was differential in any condition.
    for comp_file, comp in comparison_files.items():
        print(comp)
        comp_head = {x: i for i, x in enumerate(open(comp_file).readline().strip().split('\t'))}
        # Get a bed-object of the differential peaks to intersect with the enhancers, then get the average in case of
        # multiple overlaps. Works for both versions from DiffBind, with and w/o recentering on the summits.
        comp_fill = {x: [] for x in fill_dict}
        if not same_peaks:
            comp_peaks = []
            for entry in open(comp_file).readlines()[1:]:
                entry = entry.strip().split('\t')
                comp_peaks.append('\t'.join(entry))
            comp_peaks_inter = fill_bed.intersect(BedTool('\n'.join(comp_peaks), from_string=True), wo=True)
            # Collect the fetch_col of all peaks that intersect the enhancers.
            for inter in comp_peaks_inter:
                comp_fill['\t'.join(inter.fields[:3])] += [
                    float(inter.fields[len(fill_bed[0].fields) + comp_head[fetch_col]])]
        else:
            for entry in open(comp_file).readlines()[1:]:
                entry = entry.strip().split('\t')
                comp_fill['\t'.join(entry[:3])] = float(entry[comp_head[fetch_col]])
        # And now take the average and without a hit fill with nan. If it's the same peaks we only have one value
        for hit in comp_fill:
            fill_dict[hit][comp] = np.mean(comp_fill[hit]) if comp_fill[hit] else np.nan

    return fill_dict, list(comparison_files.values())


def promoter_fetch_col(pattern, gtf_file, tss_type='all', extend=200, gene_set=(), fetch_col='log2FC'):
    """
    Gets the promoter regions of genes and intersects them with peak files defined with a filesystem pattern to get
    their fetch_col values in the promoters. Useful for example, if one has
    multiple differential peak calls and wants to know the log2FC in promoters. If multiple regions overlap a promoter,
    their average signal is taken. Entries of promoters without overlap will be filled with NaN.

    Args:
        pattern: File path pattern e.g., DiffPeaks/DiffBind_*.bed, or the path to just one individual file.
            All files matching that pattern will be used for the intersection. The string at the asterisk will be
            used as identifier. E.g., DiffPeaks/DiffBind_Macrophages.bed will be identified as Macrophages.
            The files need to have a header starting with # with fetch_col as column name. If it's not a path pattern
            but an individual file, the suffix after the last '/' is used as identifier.
        gtf_file: gtf-file in GENCODE's format, can be gzipped.
        tss_type: "5" to get only the 5' TSS or "all" to get all unique TSS of all transcripts in the gtf-file. When 'all',
            values are averaged across multiple promoters.
        extend: Number of base pairs to extend the TSS in each direction. 200 means a window of size 401.
        gene_set: Set of Ensembl IDs or gene names or mix of both to limit the output to. If empty, return for all
            genes in the annotation.
        fetch_col: The column name in the files found by pattern to identify from which column the value should be taken.

    Returns:
        tuple:
            - **gene_values**: Nested dictionary of the Ensembl IDs as keys and the average of the fetch_col in each of the matched files.
            - **matched_files**: List of the identifiers derived from the files matching the pattern.
        """
    promoter_bed = GTF_Processing.gene_window_bed(gtf_file=gtf_file, extend=extend, tss_type=tss_type, merge=True,
                                                  gene_set=gene_set)
    prom_gene_map = {'\t'.join(x.fields[:3]): x.fields[3] for x in promoter_bed}
    fill_dict, fill_cols = peaks_fetch_col(promoter_bed, pattern, same_peaks=False, fetch_col=fetch_col)
    gene_values = {g: {c: [] for c in fill_cols} for g in set([x.fields[3] for x in promoter_bed])}
    for prom, val in fill_dict.items():
        for f_col in fill_cols:
            hit_val = float(val[f_col])
            if not np.isnan(hit_val):  # If not all promoter had a value we can form the mean after excluding NaNs.
                gene_values[prom_gene_map[prom]][f_col].append(hit_val)
    gene_values = {g: {c: np.mean(val[c]) if val[c] else np.nan for c in fill_cols} for g, val in gene_values.items()}
    return gene_values, fill_cols


def peaks_genebody_overlap(peak_file, gtf_file, gene_set=()):
    """
    Based on a bed-file path or BedTools object returns a dictionary with {gene: fraction of gene body overlapping with peak_file}.

    Args:
        peak_file: BedTool's object or path to a bed file with the regions on which the intersection is centred on.
        gtf_file: gtf-file in GENCODE's format, can be gzipped.
        gene_set: Set of Ensembl IDs or gene names or mix of both to limit the output to. If empty, return for all
            genes in the annotation.
   """
    genebody_bed = GTF_Processing.gene_body_bed(gtf_file=gtf_file, gene_set=gene_set)
    genebody_overlap = {x.fields[3]: 0 for x in genebody_bed}
    genebody_lengths = {x.fields[3]: x.length for x in genebody_bed}
    genebody_inter = genebody_bed.intersect(BedTool(peak_file).sort().merge(), wo=True)
    for inter in genebody_inter:
        genebody_overlap[inter.fields[3]] += int(inter.fields[-1])
    gene_dict = {g: genebody_overlap[g] / genebody_lengths[g] for g in genebody_overlap}
    return gene_dict

