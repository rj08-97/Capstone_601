# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')
    col = []
    marker = []
    for i in range(0, len(fold)):
        if fold[i] > 0 and prot[i] > 0:
            col.append('blue')
            marker.append('^')
        elif fold[i] < 0 and prot[i] > 0:
            col.append('orange')
            marker.append('*')
        elif fold[i] > 0 and prot[i] < 0:
            col.append('black')
            marker.append('+')
        else:
            col.append('red')
            marker.append('o')

    for i in range(len(fold)):
        # plotting the corresponding x with y
        # and respective color
        plt.scatter(fold[i], prot[i], c=col[i], s=20, marker=marker[i])
    plt.xlabel("mRNA")
    plt.ylabel("Protein")

    plt.show()

    plt.scatter(x=df["mRNA"], y=df["Protein"], c={"x": "blue", "y": "red"}, marker=["+", "*"])

    upregulated_mRNA_prot = df[(df['mRNA'] < -3) & (df['Protein'] < -5) & (df['mRNA'] > 0) & (df['Protein'] > 2)]
    upregulated_mRNA_prot.to_excel("/Users/rakshandajha/PycharmProjects/Capstone/downregulated_genes.xlsx")
    down_mRNA_prot = df[(df['mRNA'] < 0) & (df['Protein'] < -4) & (df['mRNA'] > 6) & (df['Protein'] < 3)]
    upregulated_mRNA_prot.to_excel("/Users/rakshandajha/PycharmProjects/Capstone/upregulated_genes.xlsx")
    merged_df = pd.merge(down_mRNA_prot, upregulated_mRNA_prot, on=['Gene'], how='inner')

    plt.scatter(x=down_mRNA_prot['mRNA'], y=down_mRNA_prot['Protein'], c=(down_mRNA_prot['mRNA_pvalue']))
    plt.xlabel("mRNA_downregulated")
    plt.ylabel("Protein_downregulated")
    plt.scatter(x=upregulated_mRNA_prot['mRNA'], y=upregulated_mRNA_prot['Protein'],
                c=(upregulated_mRNA_prot['mRNA_pvalue']))
    plt.xlabel("mRNA_upregulated")
    plt.ylabel("Protein_upregulated")
    plt.annotate(df['Gene'], (df['mRNA'], df['Protein']))
    plt.scatter(x=newdf['mRNA_pvalue'], y=newdf['Protein_pvalue'])

    df.plot(x="Gene", y=["mRNA_pvalue", "Protein_pvalue"], xlabel="Genes", ylabel="pvalue", legend=True)
    colors = np.random.rand(50)
    plt.scatter(x=merged_stuff["Gene"], y=merged_stuff["FoldChange_x"])
    plt.show()

    print(over_lap)
    # Replace missing values with empty strings
    # Press ⌘F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
