import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import matplotlib_venn
from matplotlib.colors import to_hex
import seaborn as sns
import scipy.stats
import itertools
import ColoursAndShapes
from Various import sanitize_path


def basic_venn(input_sets, plot_path, blob_colours=ColoursAndShapes.tol_highcontrast, title='',
               scaled=True, linestyle='', number_size=11, xsize=5, ysize=5, normalize_to=0.5, formats=['pdf']):
    """
    Based on a dictionary of {key: set} with either two or three entries do the intersection of sets and create a
    Venn diagram from it.

    Args:
        input_sets: {key: set} for all bubbles that should be plotted.
        scaled: If the bubble sizes should be scaled to the set sizes. Choose False if the size difference is too high.
        linestyle: Linestyle for the rim of the bubbles.
        normalize_to: Factor to scale the bubble size. If bubbles get cut, try decreasing it.
     """
    if 'glasbey' in blob_colours:
        blob_colours = ColoursAndShapes.glasbey_palettes[blob_colours][:len(input_sets)]
    if len(input_sets) > 3:
        print("ERROR, only making Venns for three or two sets")
        return
    key_order = list(input_sets.keys())
    if len(key_order) == 2:
        a, b = [input_sets[x] for x in key_order]
        plot_list = [len(a - b), len(b - a), len(a & b)]
    elif len(key_order) == 3:
        a, b, c = [input_sets[x] for x in key_order]
        plot_list = [len(a - b - c), len(b - a - c), len(a & b - c), len(c-a-b), len(a & c - b), len(b & c - a), len(a & b & c)]
    f, ax = plt.subplots(figsize=(xsize, ysize))
    if scaled and len(plot_list) == 3:
        v = matplotlib_venn.venn2(subsets=plot_list, set_labels=key_order, ax=ax, set_colors=blob_colours,
                                  normalize_to=normalize_to)
    elif not scaled and len(plot_list) == 3:
        v = matplotlib_venn.venn2_unweighted(subsets=plot_list, set_labels=key_order, ax=ax,
                                             set_colors=blob_colours, normalize_to=normalize_to)
    elif scaled and len(plot_list) == 7:
        v = matplotlib_venn.venn3(subsets=plot_list, set_labels=key_order, ax=ax, set_colors=blob_colours,
                                  normalize_to=normalize_to)
    elif not scaled and len(plot_list) == 7:
        v = matplotlib_venn.venn3_unweighted(subsets=plot_list, set_labels=key_order, ax=ax, set_colors=blob_colours,
                                  normalize_to=normalize_to)
    if len(plot_list) == 3:
        matplotlib_venn.venn2_circles(subsets=plot_list, linestyle=linestyle, linewidth=1, color="grey", normalize_to=normalize_to)
    elif len(plot_list) == 7:
        matplotlib_venn.venn3_circles(subsets=plot_list, linestyle=linestyle, linewidth=1, color="grey", normalize_to=normalize_to)
    for text in v.set_labels:
        if text:
            text.set_fontsize(14)
    for text in v.subset_labels:
        if text:
            text.set_fontsize(number_size)
    plt.title(title, size=16, y=1.15)
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        plt.savefig(sanitize_path(plot_path + title + '_Venn.'+form), bbox_inches='tight', format=form)
    plt.close('All')


def basic_pie(plot_df, title='', palette=None, numerate=True, legend_perc=True, output_path='', legend_title='',
              formats=['pdf']):
    """
    A pie chart.

    Args:
        plot_df: Either a DataFrame with each category as index or a dictionary with {category: count}.
        numerate: Whether to add the summed count across categories to the title.
        legend_perc: Whether the legend should write the percentages or absolute numbers.
    """
    if type(palette) == list:
        clrs = palette
    elif palette is not None and 'glasbey' not in palette:
        clrs = sns.color_palette(palette, n_colors=len(plot_df.index))
    else:
        clrs = ColoursAndShapes.glasbey_palettes[palette]
    if type(plot_df) != pd.core.frame.DataFrame:
        plot_df = pd.DataFrame.from_dict(plot_df, orient='index')
    count_list = plot_df[plot_df.columns[0]].values
    f, ax = plt.subplots()
    wedges, texts, autotexts = ax.pie(count_list, autopct="", shadow=False, startangle=0, colors=clrs,
                                      textprops={'fontsize': 11})
    if legend_perc:
        legend_labels = [str(round(count_list[i] / plot_df.sum().values[0] * 100, 2)) + '% ' + i_label for i, i_label in
                         enumerate(plot_df.index)]
    else:
        legend_labels = plot_df.index
    ax.legend(wedges, legend_labels,
              title=legend_title,
              loc="center left",
              bbox_to_anchor=(1, 0, 0.5, 1))
    ax.set_title(title + ('\n#' + str(plot_df.sum().values[0])) * numerate, fontsize=14, y=1)
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        f.savefig(sanitize_path(output_path + str(plot_df.columns[0]) + "_PieChart."+form), bbox_inches='tight',
                  format=form)
    plt.close()


def ks_test(curr_df, grouping_col, value_col):
    """
    Runs a Kolmogorov smirnov test for each non-redundant pairwise comparison. Tailored to the cumulative distribution
    function.
    """
    samples = sorted(curr_df[grouping_col].drop_duplicates().tolist())
    order_dict = {x: i for i, x in enumerate(samples)}
    comparisons = list(itertools.combinations(samples, 2))
    p_list = []
    for comp in comparisons:
        sample1_list = curr_df[curr_df[grouping_col] == comp[0]][value_col].values.tolist()
        sample2_list = curr_df[curr_df[grouping_col] == comp[1]][value_col].values.tolist()
        p_list.append([comp, scipy.stats.ks_2samp(sample1_list, sample2_list)[1]])
    p_nested = [[0 for x in range(y + 1)] for y in reversed(list(range(len(samples) - 1)))]
    p_background = [['white' for x in range(y + 1)] for y in reversed(list(range(len(samples) - 1)))]
    print(p_list)
    for i, p in enumerate(p_list):
        p_nested[order_dict[p[0][0]]][abs(order_dict[p[0][1]] - (len(samples) - 1))] = round(p[1], 5)
        if p[1] <= 0.05:
            p_background[order_dict[p[0][0]]][abs(order_dict[p[0][1]] - (len(samples) - 1))] = '#ed5739'
    for i, p in enumerate(p_nested):
        p_nested[i] += ['' for x in range((len(samples) - 1) - len(p_nested[i]))]
        p_background[i] += ['white' for x in range((len(samples) - 1) - len(p_background[i]))]

    return p_nested, p_background, samples


def cumulative_plot(plot_df, x_col, hue_col, hue_order=None, output_path='', numerate=False, title=None, xlimit=None,
                    vertical_line=False, add_all=False, table_width=0.3, table_x_pos=1.2, palette=None, font_s=16,
                    grid=True, formats=['pdf']):
    """
    Plot the cumulative distribution of all sets of grouping names in the plotting df.
    Adds a table with a Kolmogorov-Smirnov test for each non-redunant pairwise comparison next to the plot. Cells
    with a p-value ≤ 0.05 will be coloured in red.

    Args:
        plot_df: Pandas DataFrame holding the data.
        x_col: Column name in plot_df which should be plotted on the x-axis.
        hue_col: Column name in the DataFrame used for different curves.
        hue_order: If the groups in hue_col should have a specific order.
        numerate: If True, show the number of elements in each hue_col group in parentheses.
        vertical_line: To plot a vertical line at the given position, e.g. 0.
        add_all: If True, add the whole stack of values in x_col again as separate distribution labelled 'All'.
        table_width: Width of the KS table.
        table_x_pos: X-position of the KS table, to avoid it overlapping the plot.
        """
    p_df = plot_df.copy()  # Otherwise, updates for the enumeration would be referenced back to the function caller.
    if add_all:
        all_copy = p_df.copy()
        all_copy[hue_col] = 'All'
        p_df = pd.concat([all_copy, p_df], ignore_index=True)

    sns.set_style('ticks')
    sns.set_context('talk', font_scale=1.15)
    if hue_order:
        gene_set_cols = ['All']*add_all + hue_order
    else:
        gene_set_cols = p_df[hue_col].drop_duplicates().to_list()
    gene_set_cols_num = []
    different_groups = len(gene_set_cols)
    if different_groups == 3:
        linestyles = ['solid', 'dashed', 'dashdot', 'solid']
    else:
        linestyles = ['solid', 'dashed', 'solid', 'dashdot']

    if numerate:
        hue_order = []
        for i, x in enumerate(gene_set_cols):
            gene_set_cols_num.append(str(gene_set_cols[i]) + ' (n=' + str(sum(p_df[hue_col] == x)) + ')')
            hue_order.append(str(gene_set_cols[i]) + ' (n=' + str(sum(p_df[hue_col] == x)) + ')')
        for i, x in enumerate(gene_set_cols):
            p_df.loc[p_df[hue_col] == x, hue_col] = gene_set_cols_num[i]
    if not palette:
        if different_groups == 2:
            palette = ColoursAndShapes.two_contrasts[0]
        elif different_groups <= len(ColoursAndShapes.tol_vibrant):
            palette = ColoursAndShapes.tol_vibrant[:different_groups]
        else:
            palette = [to_hex(c) for c in ColoursAndShapes.glasbey_palettes['glasbey']][:different_groups]
    else:
        if 'glasbey' in palette:
            palette = [to_hex(c) for c in ColoursAndShapes.glasbey_palettes[palette]][:different_groups]
    cumu = sns.displot(p_df, x=x_col, hue=hue_col, kind='ecdf', hue_order=hue_order,
                       palette=palette, height=8, aspect=1.3, zorder=12)
    cumu.fig.subplots_adjust(top=0.93)
    cumu.fig.suptitle('Cumulative distribution of ' + x_col if not title else title, size=font_s + 6, x=0.28, y=1.01)
    cumu.ax.set_xlabel(x_col, size=font_s+4)
    cumu.ax.set_ylabel('Proportion', size=font_s+4)
    if xlimit:
        cumu.set(xlim=(xlimit[0], xlimit[1]))
    colours_order = []
    for a in range(len(cumu.ax.lines)):
        colours_order.append(cumu.ax.get_lines()[a].get_color().lower())
        cumu.ax.lines[a].set_linestyle(linestyles[a % len(linestyles)])
        cumu.ax.lines[a].set_linewidth('4')
    colour_dict = {x: i for i, x in enumerate(colours_order)}
    print(colour_dict)
    print(different_groups)
    for a in range(len(cumu.ax.lines)):
        if to_hex(cumu.legend.get_lines()[a].get_color()) in colour_dict:  # Legend has it even if there was no data.
            cumu.legend.get_lines()[a].set_linestyle(linestyles[colour_dict[to_hex(cumu.legend.get_lines()[a].get_color())] % len(linestyles)])
            cumu.legend.get_lines()[a].set_linewidth('3')
    for spine in cumu.ax.spines.values():  # Increase spline width.
        spine.set_linewidth(3)
    p_nested, p_background, samples = ks_test(p_df, hue_col, x_col)
    samples = [(str(x).split('(n=')[0][:int(len(str(x).split('(#')[0]) / 2)] + '\n' + str(x).split('(n=')[0][int(len(str(x).split('(#')[0]) / 2):]) if len(str(x).split('(n=')[0]) > 8 else x.split('(n=')[0]
               for x in samples]
    if len(samples) > 1:
        ks_table = plt.table(cellText=np.asarray(p_nested).T, cellColours=np.asarray(p_background).T,
                             rowLabels=list(reversed(samples))[:-1], colLabels=samples[:-1], loc='bottom right',
                             cellLoc='center', rowLoc='center', zorder=12, bbox=[table_x_pos, -0.11, table_width, 0.35])
        ks_table.auto_set_font_size(False)
        ks_table.set_fontsize(font_s-6)

    if grid:
        cumu.ax.grid(True, axis='both', color='#c7c7c7', linewidth=1, which='major')
    cumu.ax.set_facecolor('white')
    if vertical_line is not None:
        cumu.ax.axvline(x=vertical_line, color='#6e6e6e', zorder=1)
    if type(formats) != list:
        formats = [formats]
    for form in formats:
        cumu.savefig(sanitize_path(output_path + x_col.replace(' ', '') + '_' + hue_col + '_CumulativeDistribution.' + form), bbox_inches='tight',
                     format=form)
    plt.close()

