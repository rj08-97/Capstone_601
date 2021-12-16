# Exploring the relationship between mRNA and protein expression

Bioinformatics Advisor - Dr. Gary Skuse -> grssbi@rit.edu

Data Science Advisor - Dr. Travis Desell -> tjdvse@rit.edu


The goals of this study are to 1) explore relationships between differentially expressed mRNAs and their encoded proteins to determine whether there is a linear relationship between mRNA expression levels and the levels of their encoded proteins, 2) identify relationships among genes carrying mutations and the pathways in which they play a role to determine the effect those mutations have on complex biological pathways and 3) identify relationships among the expression levels of genes that encode transcription factors and the genes they regulate. The important questions that will be addressed in this study are, 1) does every transcribed mRNA leads to the translation of protein, 2) how do the identified mutations in transcription factors, gene products that regulate transcription, affect the process of transcription and translation. Moreover, the study will also focus on the relationship between the transcription factors and the mutations in these factors will be studied. For example, if there is a rise in the mutation, does that affect the transcription process. 


Biology is complicated, so much that happens within the human body that we do not  understand. There was a time when biology was considered simple science, and it was not until the discovery of DNA structure that made scientists realize just how much we do not know. Gene expression and regulation describes the complexity of a human body, transcription is an important stage in gene expression. Therefore, having a good understanding of the relationship between transcription and translation is very important as they play a very significant role in every cell within the human body.  Disruption of either process can lead to disease. 


# How to run the code

The code is pretty straight forward and well commented the only package required is the biopython to run the scripts above. The first step is to run the preliminary analysis to perform and visualize the correlation results. Then is to run the data mining file to gather all the information from the genomic sequences like sequence alignment, similarity score, sequence complexity, motif analysis and transcription factors. Lastly, is to run the RNA secondary structures file in the folder of mRNA folds to generate the dot structures and gain the free energy of each sequences. Moreover, it produces the RNA binding sites and its mapping within the secondary structure.

# Dataset

The was from the study conducted by Fernandez4 et al. [1] where they use and express the mRNA and protein expression from 14 cancer cell lines from Low-grade Serous Ovarian Carcinoma. The dataset is included with the paper and can be downloaded from the link below [1]. Moreover, the genomic sequences are retrived from the NCBI database.


# Citations

[1] Shrestha R, Llaurado Fernandez M, Dawson A, Hoenisch J, Volik S, Lin YY, Anderson S, Kim H, Haegert AM, Colborne S, Wong NKY, McConeghy B, Bell RH, Brahmbhatt S, Lee CH, DiMattia GE, Le Bihan S, Morin GB, Collins CC, Carey MS. Multiomics Characterization of Low-Grade Serous Ovarian Carcinoma Identifies Potential Biomarkers of MEK Inhibitor Sensitivity and Therapeutic Vulnerability. Cancer Res. 2021 Apr 1;81(7):1681-1694. doi: 10.1158/0008-5472.CAN-20-2222. Epub 2021 Jan 13.
PMID: 33441310.
https://pubmed.ncbi.nlm.nih.gov/33441310/
