import pandas as pd
import matplotlib as mat
import matplotlib.pyplot as plt
import numpy as np
import seaborn
from scipy.stats import pearsonr, spearmanr


def correlation(merge_df):
    """
    The method performs correlation testing on mRNA and protein expression data and plots heatmaps.
    :param merge_df: The overlap genes of mRNA and protein
    :return: correlation coefficients and heat maps
    """

    df = merge_df.rename(columns={'FoldChange_x': 'mRNA', 'FoldChange_y': 'Protein','qvalue_x' : 'mRNA_qvalue',
                                      'pvalue_x': 'mRNA_pvalue', 'pvalue_y': 'Protein_pvalue', 'qvalue_y' : 'Protein_qvalue',})

    # To find the correlation among
    # the columns using pearson method
    correlation_ = df.corr(method='pearson')

    # the columns using kendall method
    correlation_kendal = df.corr(method='kendall')

    # the columns using spearman method
    correlation_spear = df.corr(method='spearman')

    protein = pd.DataFrame({"Gene": df["Gene"],"Protein": df["Protein"]})
    mRNA = pd.DataFrame({"Gene": df["Gene"],"mRNA": df["mRNA"]})

    corr_df = mRNA.corrwith(protein, method='pearson')
    print(1)
    seaborn.heatmap(correlation_)
    plt.show()

    # calculate Pearson's correlation
    corr, _ = pearsonr(df["mRNA"], df["Protein"])
    print('Pearsons correlation: %.3f' % corr)

    # calculate spearman's correlation
    corr_spearmen, _ = spearmanr(df["mRNA"], df["Protein"])
    print('Spearmans correlation: %.3f' % corr_spearmen)
    print(1)

def overlap(mRNA, protein):
    """
    Finding the genes overlap between the mRNA and protein dataset and performing preliminary data visualization
    :param mRNA: The mRNA expression data with fold change and p-value
    :param protein: The protein expression data with fold change and p-value
    :return: The overlap genes between mRNA and protein along with prelim data visualization.
    """
    merged_stuff = pd.merge(mRNA, protein, on=['Gene'], how='inner')
    over_lap = mRNA['Gene'].isin(protein['Gene']).value_counts()

    # renaming the columns for to avoid mistake

    df = merged_stuff.rename(columns={'FoldChange_x': 'mRNA', 'FoldChange_y': 'Protein',
                                      'pvalue_x': 'mRNA_pvalue', 'pvalue_y': 'Protein_pvalue'})

    df.to_excel("/Users/rakshandajha/PycharmProjects/Capstone/overlap.xlsx")

    newdf = df[(df['mRNA'] > -5) & (df['Protein'] > -5)]
    # newdf.to_excel("/Users/rakshandajha/PycharmProjects/Capstone/overlap_genes.xlsx")

    protein = pd.DataFrame({"Gene": df["Gene"], "Protein": df["Protein"]})
    mRNA = pd.DataFrame({"Gene": df["Gene"], "mRNA": df["mRNA"]})

    fold = np.array(newdf["mRNA"])
    prot = np.array(newdf["Protein"])



if __name__ == "__main__":
    # Load data
    mRNA = pd.read_excel("/Users/rakshandajha/PycharmProjects/Capstone/mRNA.xlsx")
    protein = pd.read_excel("/Users/rakshandajha/PycharmProjects/Capstone/Protein.xlsx")
    merged_stuff = pd.merge(mRNA, protein, on=['Gene'], how='inner')
    over_lap = mRNA['Gene'].isin(protein['Gene']).value_counts()
    correlation(merged_stuff)
    overlap(mRNA, protein)
