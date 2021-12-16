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


def readingFile():
    """The function reads in the sequences in fasta format and returns it as a dictionary"""
    seq = SeqIO.parse(open("/Users/rakshandajha/Documents/Capstone_601/DNA_Seq/GRCh38_latest_rna.fasta"), 'fasta')
    fasta_name = []
    fasta_sequence = []
    for s in seq:
         fasta_name.append(s.id)
         fasta_sequence.append(str(s.seq))
    data_seq = []
    for i in range(0, len(seq)):
         if seq[i] in fasta_name:
            data_seq.append(fasta_sequence[i])

    return data_seq


def pairwise(sequence):
    """The function produces the local and global alignment along with the similarity scores."""
    sequence = list(sequence.values())
    align = []
    score = []
    l_score = []
    local_align = []
    for i in range(0, len(sequence) - 1):
        for j in range(0, len(sequence)):
            align = pairwise2.align.globalxx(sequence[i], sequence[j])
            seq_length = min(len(sequence[i]), len(sequence[j]))
            matches = align[0][2]
            score = (matches / seq_length) * 100
            # print("The sequence similarity score for global alignment ", score)
            local_align = pairwise2.align.localxx(sequence[i], sequence[j])
            local_matches = local_align[0][2]
            l_score.append((local_matches / seq_length) * 100)

    print("Needleman-Wunsch Global Algorithm")
    for index in range(len(align)):
        print(format_alignment(*align[index]))
    print("Smith-Waterman Local Algorithm")
    for index in range(len(local_align)):
        print(format_alignment(*local_align[index]))
        print("The sequence similarity score for local alignment ", l_score[index])
    return np.array(align), np.array(local_align), score, l_score


def counts(sequence, start, end):
    """The function produces the count and position of start and stop codons."""
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
    # sns.distplot(count_stop_sequences, bins="doane", axlabel="Frequency of start and stop codon", kde=False,  hist_kws={"align": "right"})
    # fig.legend(labels=['Start codons', 'Stop codons'])
    # plt.show()

    return count_start_sequence, count_stop_sequences


def motifAnalysis(sequence, start):
    """The function produces the motifs and its frequency along with its position"""
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
    """The following funtion encodes the data in one hot encoding"""
    input = list(input.values())
    sequence = []
    mapping = dict(zip("ACGT", range(4)))
    for index in range(len(input)):
        sequence = [mapping[i] for i in input[index]]

    sequence = np.eye(4)[sequence]
    print(sequence)
    return sequence


def validation(sequence):
    """The following function validates the input sequences"""
    dna_sequences = list(sequence.values())
    dna_char = "ACGT"
    for i in range(len(dna_sequences)):
        for letter in dna_sequences[i]:
            if letter not in dna_char:
                print("Not a valid sequence")
                return False
        return True


def GCContent(data):
    """The following function performs a count on the number of G and C present in the sequence"""
    data = list(data.values())
    GC = []
    gc = []
    for i in range(0, len(data)):
        sequence = Seq(data[i])
        GC.append(round((sequence.count('G') + sequence.count('C')) / len(sequence) * 100))
        gc.append(SeqUtils.GC_skew(sequence, 100))

    print(gc)
    return GC


def main():
    # mRNA = pd.read_excel("/Users/rakshandajha/Documents/Capstone_601/DNA_Seq/overlap.xlsx")
    start = "ATG"
    end = {1: "TAG", 2: "TAA", 3: "TGA"}
    # getSequence(mRNA)
    input_sequences = readingFile()
    print("validation step")
    validation(input_sequences)
    print("GC content")
    GCContent(input_sequences)
    count_start_sequence, count_stop_sequences = counts(input_sequences, start, end)
    print(count_start_sequence, ' ',  count_stop_sequences)
    print("Alignment")
    alignment = pairwise(input_sequences)
    print("Motif Analysis")
    start, weight = motifAnalysis(input_sequences, start)

main()
