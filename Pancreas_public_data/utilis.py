import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
from datetime import date

today = date.today()

def helper_violin_plot(data, colors):
    '''
    This function creates a violin plot with a boxplot inside.
    
    Parameters
    ----------
    data : list
        List of lists with data to plot.
    colors : list
        List of colors to use for the plot.
    '''
    
    plt.figure(figsize = (2, 2.5), dpi = 300)
    
    colors = colors
    customPalette = sns.set_palette(sns.color_palette(colors))
    
    ax = sns.violinplot(data = data, saturation = 0.9, width = 0.9, palette = customPalette, linewidth = 0.3, kws = {'linecolor' : 'black'})
    for i, c in enumerate(ax.collections):
        ax.collections[i].set_edgecolor('black')

    sns.boxplot(data = data, width = 0.4,
                boxprops = {'zorder': 2, 'edgecolor' : 'black'},
                capprops = {'color' : 'black'},
                whiskerprops = {'color' : 'black'},
                medianprops = {'color' : 'black'},
                showfliers = False,
                linewidth = 0.3,
                ax = ax)

                                                        
    sns.stripplot(data = data, color = 'black', ax = ax, size = 0.4)
    return ax



def violin_plot(original_adata, colors, group, group_conditions, alternatives = ['two-sided', 'greater', 'less'], umap_color = 'Purples', marker = None, adata_gene = None, tissue = None, geneset = None, geneset_name = None, module_score = False, n_bins = 25, ctrl_size = 50):

    data = {}

    if module_score:
        score_name = f'module score for {geneset_name}, {group_conditions}'
        sc.tl.score_genes(original_adata, geneset, score_name = score_name, ctrl_size = ctrl_size, n_bins = n_bins)

        data[group_conditions[0]] = np.array(original_adata.obs[original_adata.obs[group] == group_conditions[0]][score_name])
        data[group_conditions[1]] = np.array(original_adata.obs[original_adata.obs[group] == group_conditions[1]][score_name])
        
        y_label = 'Module score'
        title = f'Module Score for {" vs. ".join(group_conditions)}\nbased on geneset - {geneset_name}'
        umap_title = f'Module Score based on geneset - {geneset_name}'
        umap_group = score_name
    else:
        adata_0 = adata_gene[adata_gene.obs[group] == group_conditions[0]]
        adata_1 = adata_gene[adata_gene.obs[group] == group_conditions[1]]
        
        data[group_conditions[0]] = list(adata_0.X[:, np.where(adata_0.var_names == marker)])
        data[group_conditions[1]] = list(adata_1.X[:, np.where(adata_1.var_names == marker)])

        y_label = 'Log-scaled expression'
        title = f'Log-scaled expression of {marker} in {group}'
        umap_title = f'Log-scaled expression of {marker}'
        umap_group = marker
    
    p_values = []
    for alternative in alternatives:
        _ , p = stats.ranksums(data[group_conditions[0]], data[group_conditions[1]], alternative = alternative)
        p_values.append(p)
    p_values = [f'{i[0]} : {round(i[1], 6)}' for i in list(zip(alternatives, p_values))]
    
    ax = helper_violin_plot(list(data.values()), colors)

    ax.set_yticklabels(ax.get_yticks(), size = 4)    
    ax.set_xticklabels(group_conditions, size = 4)
    
    ax.set_title(f'{title}\n\nwilcoxon rank sum, alternative hypothesis : p value\n{", ".join(p_values)}', fontsize = 4)
    sns.despine()

    plt.savefig('Pancreas_EX0008_violin_plot.pdf', dpi = 300, bbox_inches = 'tight')
    plt.show()
    plt.clf()

    #sc.pl.umap(original_adata, color = umap_group, title = umap_title, show = False, color_map = umap_color)    
    #plt.savefig(f'figures/{today}/{today}_MAA_EX0008_umap_{umap_title}.pdf', dpi = 300, bbox_inches = 'tight')
    #plt.show()
    #plt.clf()
        
 