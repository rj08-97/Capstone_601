from Bio import SeqIO
import pandas as pd
import numpy as np
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqUtils
from Bio import motifs
from Bio import Entrez
from pydna.genbank import Genbank
from pyensembl import EnsemblRelease
import pyseqlab


def readingFile(data):
    global avil_sequence, anxa1_sequence, bst2_sequence, cpne4_sequence, cryab_sequence, \
        gpd2_sequence, hmga2_sequence, krt16_sequence, krt17_sequence, \
        lcn2_sequence, mtap_sequence, epcam_sequence, avil_name, bst2_name, cpne4_name, cryab_name, epcam_name, hmga2_name, gpd2_name, krt16_name, lcn2_name, krt17_name, mtap_name
    anxa1_sequences = SeqIO.parse(open("../DNA_Seq/ANXA1_mRNA.fasta"), 'fasta')
    avil_sequences = SeqIO.parse(open("../DNA_Seq/AVIL_mRNA.fasta"), 'fasta')
    bst2_sequences = SeqIO.parse(open("../DNA_Seq/BST2_mRNA.fasta"), 'fasta')
    cpne4_sequences = SeqIO.parse(open("../DNA_Seq/CPNE4_mRNA.fasta"), 'fasta')
    cryab_sequences = SeqIO.parse(open("../DNA_Seq/CRYAB_mRNA.fasta"), 'fasta')
    epcam_sequences = SeqIO.parse(open("../DNA_Seq/EPCAM_mRNA.fasta"), 'fasta')
    gpd2_sequences = SeqIO.parse(open("../DNA_Seq/GPD2_mRNA.fasta"), 'fasta')
    hmga2_sequences = SeqIO.parse(open("../DNA_Seq/HMGA2_mRNA.fasta"), 'fasta')
    krt16_sequences = SeqIO.parse(open("../DNA_Seq/KRT16_mRNA.fasta"), 'fasta')
    krt17_sequences = SeqIO.parse(open("../DNA_Seq/KRT17_mRNA.fasta"), 'fasta')
    lcn2_sequences = SeqIO.parse(open("../DNA_Seq/LCN2_mRNA.fasta"), 'fasta')
    mtap_sequences = SeqIO.parse(open("../DNA_Seq/MTAP_mRNA.fasta"), 'fasta')

    for fasta in anxa1_sequences:
        anxa1_name, anxa1_sequence = fasta.id, str(fasta.seq)
        # new_sequence = some_function(sequence)
    for sequ in avil_sequences:
        avil_name, avil_sequence = sequ.id, str(sequ.seq)
        # new_sequence = some_function(sequence)
    for fasta in bst2_sequences:
        bst2_name, bst2_sequence = fasta.id, str(fasta.seq)
        # new_sequence = some_function(sequence)
    for fast in cpne4_sequences:
        cpne4_name, cpne4_sequence = fast.id, str(fast.seq)
        # new_sequence = some_function(sequence)
    for fasta in cryab_sequences:
        cryab_name, cryab_sequence = fasta.id, str(fasta.seq)
        # new_sequence = some_function(sequence)
    for fasta in epcam_sequences:
        epcam_name, epcam_sequence = fasta.id, str(fasta.seq)
        # new_sequence = some_function(sequence)
    for fasta in gpd2_sequences:
        gpd2_name, gpd2_sequence = fasta.id, str(fasta.seq)
        # new_sequence = some_function(sequence)
    for fasta in hmga2_sequences:
        hmga2_name, hmga2_sequence = fasta.id, str(fasta.seq)
        # new_sequence = some_function(sequence)
    for fasta in krt16_sequences:
        krt16_name, krt16_sequence = fasta.id, str(fasta.seq)
        # new_sequence = some_function(sequence)
    for fasta in krt17_sequences:
        krt17_name, krt17_sequence = fasta.id, str(fasta.seq)
        # new_sequence = some_function(sequence)
    for fasta in lcn2_sequences:
        lcn2_name, lcn2_sequence = fasta.id, str(fasta.seq)
        # new_sequence = some_function(sequence)
    for fasta in mtap_sequences:
        mtap_name, mtap_sequence = fasta.id, str(fasta.seq)
        # new_sequence = some_function(sequence)

    input_name = [avil_name, avil_name, bst2_name, cpne4_name, cryab_name, epcam_name, gpd2_name, hmga2_name,
                  krt16_name, krt17_name,
                  lcn2_name, mtap_name]
    seq = SeqIO.parse(open("/Users/rakshandajha/Documents/Capstone_601/DNA_Seq/GRCh38_latest_rna.fasta"), 'fasta')
    fasta_name = []
    fasta_sequence = []
    for s in seq:
        fasta_name.append(s.id)
        fasta_sequence.append(str(s.seq))
    data_seq = []
    for i in range(0, len(input_name)):
        if input_name[i] in fasta_name:
            data_seq.append(fasta_sequence[i])

    input_sequences = {"ANXA1": anxa1_sequence, "AVIL": avil_sequence,
                       "BST2": bst2_sequence, "CPNE4": cpne4_sequence, "CRYAB": cpne4_sequence,
                       "EPCAM": epcam_sequence, "GPD2": gpd2_sequence, "HMGA2": hmga2_sequence,
                       "KRT16": krt16_sequence, "KRT17": krt17_sequence, "LCN2": lcn2_sequence, "MTAP": mtap_sequence
                       }

    return input_sequences, data_seq


def pairwise(sequence, dna_seq):
    sequence = list(sequence.values())
    align = []
    score = []
    l_score = []
    local_align = []
    for i in range(0, len(sequence) - 1):
        align = pairwise2.align.globalxx(sequence[i], dna_seq[i])
        seq_length = min(len(sequence[i]), len(dna_seq[i]))
        matches = align[0][2]
        score = (matches / seq_length) * 100
        # print("The sequence similarity score for global alignment ", score)
        local_align = pairwise2.align.localxx(sequence[i], dna_seq[i])
        local_matches = local_align[0][2]
        l_score.append((local_matches / seq_length) * 100)

    # print("Needleman-Wunsch Global Algorithm")
    # for index in range(len(align)):
    # print(format_alignment(*align[index]))
    print("Smith-Waterman Local Algorithm")
    for index in range(len(local_align)):
        print(format_alignment(*local_align[index]))
        print("The sequence similarity score for local alignment ", l_score[index])
    return np.array(align), np.array(local_align), score, l_score


def counts(sequence, start, end):
    sequence = list(sequence.values())
    stop_seq = list(end.values())
    seq = ""
    count_start_sequence = []
    count_TAG_sequence = []
    count_TAA_sequence = []
    count_TGA_sequence = []
    count_stop_sequences = []
    for i in range(len(sequence)):
        seq = Seq(sequence[i])
        count_start_sequence.append(seq.find(start))
        count_TAG_sequence.append(seq.find(stop_seq[0]))
        count_TAA_sequence.append(seq.find(stop_seq[1]))
        count_TGA_sequence.append(seq.find(stop_seq[2]))
    count_stop_sequences = pd.DataFrame(
        {'TAG': count_TAG_sequence,
         'TAA': count_TAA_sequence,
         'TGA': count_TGA_sequence
         })
    fig = sns.distplot(count_start_sequence, bins="doane", axlabel="Frequency of start and stop codon", kde=False,
                       hist_kws={"align": "right"})
    sns.distplot(count_stop_sequences, bins="doane", axlabel="Frequency of start and stop codon", kde=False,
                 hist_kws={"align": "right"})
    fig.legend(labels=['Start codons', 'Stop codons'])
    plt.show()

    return count_start_sequence, count_stop_sequences


def motifAnalysis(sequence, start):
    pattern = Seq(start)
    start_results = SeqUtils.nt_search(str(sequence), pattern)
    instances = [Seq("TTTGTA"), Seq("AATAAA"), Seq("CCTCCC"), Seq("TGTGTG"), Seq("TTATTT"), Seq("CAGTTT")]
    motif = motifs.create(instances)
    print(motif.degenerate_consensus)
    print(motif.counts)
    motif.weblogo("giant_LOGO.pdf", format="pdf")
    weightMatrix = motif.pssm
    print(weightMatrix)
    return start_results, weightMatrix


def dataEncoding(input):
    input = list(input.values())
    sequence = []
    mapping = dict(zip("ACGT", range(4)))
    for index in range(len(input)):
        sequence = [mapping[i] for i in input[index]]

    sequence = np.eye(4)[sequence]
    print(sequence)
    return sequence


def validation(sequence):
    dna_sequences = list(sequence.values())
    dna_char = "ACGT"
    for i in range(len(dna_sequences)):
        for letter in dna_sequences[i]:
            if letter not in dna_char:
                print("Not a valid sequence")
                return False
        return True


def GCContent(data):
    data = list(data.values())
    GC = []
    gc =[]
    for i in range(0, len(data)):
        sequence = Seq(data[i])
        GC.append(round((sequence.count('G') + sequence.count('C')) / len(sequence) * 100))
        gc.append( SeqUtils.GC_skew(sequence, 100))

    print(gc)
    return GC


def main():
    mRNA = pd.read_excel("/Users/rakshandajha/Documents/Capstone_601/DNA_Seq/overlap.xlsx")
    start = "ATG"
    # getSequence(mRNA)
    input_sequences, dna_seq = readingFile(mRNA)
    validation(input_sequences)
    GCContent(input_sequences)
    # counts(input_sequences, start, end)
    alignment = pairwise(input_sequences, dna_seq)
    start, weight = motifAnalysis(input_sequences, start)
    encoded_seq = dataEncoding(input=input_sequences)


main()
