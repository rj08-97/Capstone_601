from Bio import SeqIO
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from numpy import argmax

import DNA_Seq
from sklearn.preprocessing import OneHotEncoder, LabelEncoder


def readingFile():
    global avil_sequence, anxa1_sequence, bst2_sequence, cpne4_sequence, cryab_sequence, \
        gpd2_sequence, hmga2_sequence, krt16_sequence, krt17_sequence, \
        lcn2_sequence, mtap_sequence, epcam_sequence
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

    input_sequences = {"ANXA1": anxa1_sequence, "AVIL": avil_sequence,
                       "BST2": bst2_sequence, "CPNE4": cpne4_sequence, "CRYAB": cpne4_sequence,
                       "EPCAM": epcam_sequence, "GPD2": gpd2_sequence, "HMGA2": hmga2_sequence,
                       "KRT16": krt16_sequence, "KRT17": krt17_sequence, "LCN2": lcn2_sequence, "MTAP": mtap_sequence
                       }

    return input_sequences


def counts(sequence, start,end):
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

    return count_start_sequence, count_stop_sequences


def dataEncoding(input, orf):
    input = list(input.values())
    orf_input = list(orf.values())
    sequence = []
    mapping = dict(zip("ACGT", range(4)))
    for index in range(len(input)):
        sequence = [mapping[i] for i in input[index]]
    orf_seq = []
    for index in range(len(orf_input)):
        orf_seq = [mapping[ind] for ind in orf_input[index]]
    sequence = np.eye(4)[sequence]
    orf_seq = np.eye(4)[orf_seq]
    return sequence, orf_seq


def main():
    start = "ATG"
    end = {1:"TAG",2: "TAA",3:" TGA"}
    Orf = {
        "ANXA1": "ATGGCAATGGTATCAGAATTCCTCAAGCAGGCCTGGTTTATTGAAAATGAAGAGCAGGAATATGTTCAAACTGTGAAGTCATCCAAAGGTGGTCCCGGATCAGCGGTGAGCCCCTATCCTACCTTCAATCCATCCTCGGATGTCGCTGCCTTGCATAAGGCCATAATGGTTAAAGGTGTGGATGAAGCAACCATCATTGACATTCTAACTAAGCGAAACAATGCACAGCGTCAACAGATCAAAGCAGCATATCTCCAGGAAACAGGAAAGCCCCTGGATGAAACACTGAAGAAAGCCCTTACAGGTCACCTTGAGGAGGTTGTTTTAGCTCTGCTAAAAACTCCAGCGCAATTTGATGCTGATGAACTTCGTGCTGCCATGAAGGGCCTTGGAACTGATGAAGATACTCTAATTGAGATTTTGGCATCAAGAACTAACAAAGAAATCAGAGACATTAACAGGGTCTACAGAGAGGAACTGAAGAGAGATCTGGCCAAAGACATAACCTCAGACACATCTGGAGATTTTCGGAACGCTTTGCTTTCTCTTGCTAAGGGTGACCGATCTGAGGACTTTGGTGTGAATGAAGACTTGGCTGATTCAGATGCCAGGGCCTTGTATGAAGCAGGAGAAAGG"
                 "AGAAAGGGGACAGACGTAAACGTGTTCAATACCATCCTTACCACCAGAAGCTATCCACAACTTCGCAGAGTGTTTCAGAAATACACCAAGTACAGTAAGCATGACATGAACAAAGTTCTGGACCTGGAGTTGAAAGGTGACATTGAGAAATGCCTCACAGCTATCGTGAAGTGCGCCACAAGCAAACCAGCTTTCTTTGCAGAGAAGCTTCATCAAGCCATGAAAGGTGTTGGAACTCGCCATAAGGCATTGATCAGGATTATGGTTTCCCGTTCTGAAATTGACATGAATGATATCAAAGCATTCTATCAGAAGATGTATGGTATCTCCCTTTGCCAAGCCATCCTGGATGAAACCAAAGGAGATTATGAGAAAATCCTGGTGGCTCTTTGTGGAGGAAACTAA",
        "AVIL": 'ATGCCTCTGACCAGTGCCTTCAGGGCTGTGGACAACGACCCTGGGATCATTGTCTGGAGAATAGAGAAAATGGAGCTGGCGCTGGTGCCTGTGAGCGCCCACGGCAACTTCTATGAGGGGGACTGCTACGTCATCCTCTCGACCCGGAGAGTGGCCAGTCTCCTATCCCAGGACATCCACTTCTGGATCGGGAAGGACTCCTCCCAGGATGAGCAAAGCTGCGCAGCCATATATACCACACAGCTGGACGACTACCTGGGAGGCAGCCCTGTGCAGCACCGAGAGGTCCAGTACCATGAGTCAGACACTTTCCGTGGCTACTTCAAGCAGGGCATCATCTACAAGCAGGGGGGTGTCGCCTCTGGGATGAAGCACGTGGAGACCAATACCTACGACGTGAAGCGGCTG'
                'CTACATGTGAAAGGGAAAAGAAACATCAGGGCTACCGAGGTGGAAATGAGCTGGGACAGTTTCAACCGAGGTGATGTCTTCTTGCTGGACCTTGGGAAAGTCATCATCCAATGGAATGGCCCAGAGAGCAACAGTGGGGAGCGCCTGAAGGCTATGCTTCTGGCAAAGGATATTCGAGACAGGGAGCGAGGGGGCCGTGCTAAAATAGGAGTGATCGAGGGAGACAAGGAGGCAGCCAGCCCAGAGCTGATGAAGGTCCTTCAGGACACCCTTGGCCGACGCTCCATTATCAAGCCTACAGTCCCTGATGAGATCATAGATCAGAAGCAGAAATCAACTATCATGTTGTATCATATCTCAGATTCAGCTGGGCAGCTGGCAGTCACAGAGGTAGCAACAAGGCCTCTGGTCCAGGACTTACTGAACCATGATGACTGCTACA'
                'TCCTGGACCAAAGTGGAACCAAAATCTACGTGTGGAAAGGAAAAGGAGCCACAAAGGCTGAAAAACAGGCAGCCATGTCTAAAGCGCTGGGCTTCATCAAGATGAAGAGCTACCCCAGCAGCACCAATGTGGAGACCGTCAACGATGGTGCTGAGTCGGCCATGTTCAAGCAGCTGTTCCAGAAGTGGTCAGTAAAGGACCAGACCATGGGCCTGGGGAAAACGTTCAGCATTGGTAAAATTGCTAAAGTTTTCCAGGATAAATTTGATGTGACTCTGCTACACACCAAGCCAGAGGTAGCTGCCCAGGAAAGAATGGTCGATGATGGCAACGGAAAAGTTGAGGTCTGGAGAATTGAGAACCTGGAGCTGGTCCCTGTGGAGTATCAATGGTATGGCTTCTTTTATGGGGGAGACTGTTATCTGGTCCTCTACACATACGAGGTAAATGGGAAGCCACATCACATCTTGTACATCTGGCAGGGCCGCCACGCCTCACAGGATGAGCTGGCAGCCTCAGCATACCAGGCAGTGGAGGTGGATCGGCAGTTTGATGGGGCTGCTGTGCAGGTTCGAGTCAGGATGGGAACGGAGCCACGCCACTTCATGGCCATCTTCAAAGGGAAGCTAGTTATCTTTGAGGGTGGGA'
                'CTTCCAGGAAGGGAAATGCCGAGCCTGACCCTCCAGTAAGACTCTTCCAAATTCATGGAAATGACAAATCTAACACCAAAGCAGTGGAAGTTCCAGCCTTTGCCTCCTCCCTAAACTCCAATGATGTCTTTCTGCTGCGAACTCAGGCAGAGCACTACCTGTGGTATGGCAAGGGGTCTAGTGGGGATGAGCGGGCAATGGCTAAGGAGCTGGCCAGCCTTCTCTGTGATGGCAGCGAGAACACTGTGGCCGAGGGCCAGGAGCCAGCCGAGTTCTGGGACCTACTGGGAGGGAAAACTCCCTATGCCAATGATAAAAGACTTCAGCAGGAAATCCTAGATGTCCAGTCTCGTCTCTTTGAATGTTCCAATAAGACCGGCCAATTCGTTGTCACTGAGATCACAGACTTCACCCAGGATGACCTGAACCCTACTGACGTGATGCTCCTAGATACCTGGGACCAGGTGTTCTTGTGGATTGGGGCTGAGGCCAATGCCACGGAGAAGGAGAGTGCCCTTGCCACAGCACAGCAGTACCTGCACACTCACCCCAGCGGCCGAGATCCCGACACACCAATCCTGATCATTAAGCAGGGGTTTGAGCCTCCCATCTTCACAGGCTGGTTCCTAGCCTGGGACCCTAACATTTGGAGTGCAGGAAAAACATATGAACAATTAAAAGAAGAGCTGGGAGATGCTGCTGCTATCATGCGAATCACTGCTGACATGAAGAATGCAACCCTCTCCCTGAATTCTAATGACAGTGAGCCAAAATATTACCCTATAGCAGTTCTGTTGAAAAACCAGAATCAGGAGCTGCCTGAGGATGTAAACCCTGCCAAAAAGGAGAATTACCTCTCTGAACAGGACTTTGTGTCTGTGTTTGGCATCACAAGAGGGCAATTTGCAGCTCTGCCTGGCTGGAAACAGCTCCAAATGAAGAAAGAAAAGGGGCTTTTCTAA',
        "BST2": "ATGGCATCTACTTCGTATGACTATTGCAGAGTGCCCATGGAAGACGGGGATAAGCGCTGTAAGCTTCTGCTGGGGATAGGAATTCTGGTGCTCCTGATCATCGTGATTCTGGGGGTGCCCTTGATTATCTTCACCATCAAGGCCAACAGCGAGGCCTGCCGGGACGGCCTTCGGGCAGTGATGGAGTGTCGCAATGTCACCCATCTCCTGCAACAAGAGCTGACCGAGGCCCAGAAGGGCTTTCAGGATGTGGAGGCCCAGGCCGCCACCTGCAACCACACTGTGATGGCCCTAATGGCTTCCCTGGATGCAGAGAAGGCCCAAGGACAAAAGAAAGTGGAGGAGCTTGAGGGAGAGATCACTACATTAAACCATAAGCTTCAGGACGCGTCTGCAGAGGTGGAGCGACTGAGAAGAGAAAACCAGGTCTTAAGCGTGAGAATCGCGGACAAGAAGTACTACCCCAGCTCCCAGGACTCCAGCTCCGCTGCGGCGCCCCAGCTGCTGATTGTGCTGCTGGGCCTCAGCGCTCTGCTGCAGTGA",
        "CPNE4": "ATGTTTCAAAATGCTTTAGAGTGTAAAGATAAGAAATTTAGCAAAATTTCAACAAACCTGTATAAAAAATGGATACAGATACCATTTTATTTCAAGGCGACATGCTGTAATGTAAATAGTAAGTTATTGTTATCTGTGTTTCTTTGTAAAAAGTTAACAATACAAGGGCTTAACAGATCCAGGCCTCAGGGCTACACATGCAAAAAAAATTGTCAGATAGTTAAAAATCACCAAAACGTGCTATTTTTAAATGTGTATATGTTGTTGGTTTTTTAA",
        "CRYAB": "ATGGGAATGGTGCGCTCAGGGCCAGAGACCTGTTTCCTTGGTCCATTCACAGTGAGGACCCCATCAGATGACAGGGATGAAGTAATGGTGAGAGGGTCTACATCAGCTGGGATCCGGTATTTCCTGTGGAACTCCCTGGAGATGAAACCATGTTCATCCTGGCGCTCTTCATGTTTTCCATGCACCTCAATCACATCTCCCAACACCTTAACTTTGAGTTCCTCTGGGGAGAAGTGCTTCACATCCAGGTTGACAGAGAACCTGTCCTTCTCCAGGCGCATCTCTGAGAGTCCAGTGTCAAACCAGCTGGGTGCCCGCAGGAAGGAGGGTGGCCGAAGGTAG",
        "EPCAM": "ATGCTGCGCGCGCGCCGAGAAGAGGGGCGCGGGAGGGGCCCGGGACTCGCGGGAGGGCGTGCGCCGGGAGTGGGACACGAGGAGCCGGCCGGGCAGCGCGAGGCCTGGGGCACGCGGGTCCGCGTCGGGAGGACAGCGACGAGGGGGTCCCCGGACCGCGTCGAAGGTGCTCGCTCGCCGAAGGACTAG",
        "GPD2": "ATGTCCTGGCAGCATGGAGTGGAATCCGTCCTCTTGTTACAGACCCCAAATCTGCAGATACTCAGTCTATCTCCCGAAATCATGTTGTTGATATCAGTGAGAGTGGCCTTATTACTATAG",
        "HMGA2": "ATGGGGAGACTCCGCCGGGATAGGGTCGGGCACGGAGCACAGGCAGAGGACAGAGTAGTGGGTGGCACCGCGCCTCCTAGGGTGGCGGGAGCAGGCAGCGGCGGCACCGGAGAGTCGGAGGGGGACGGGCTGCTAGCTCCTGAGTCTTGCACCAAGCGCGCGCTGCCCGGGCTGGAAGTTTTCTGA",
        "KRT16": "ATGCTTGCTGGGAGGAAAGGTGGGCATCCTCGCCCTCCAGCAGGCGGCGGTAGGTGGCAATCTCCTGCTCCAGCCGCGTCTTCACATCCAGCAAGATCTGGTACTCCTGGCTCTGCTGCTCCATCTCACAGCGTAGCTGGGCCAGCTGCTCCTCCACACTGCCAATCAGTCCCTGGATCTGGGACAGCTGCATGCAGTAGCGGCCTTTGGTCTCCTCCAGGCTGTTCTCCAGGGATGCTTTCATGCTGA",
        "KRT17": "ATGCTGAGCTGGGACTGCAGCTCTATCTCCAAGGCCTGCATGGTGCGCCGGAGCTCCGAGATCTCACTCTTGCCACTCTGCACCAGCTCACTGTTGGTGGCCACCTCGCGGTTCAGTTCCTCTGTCTTGCTGAAGAACCAATCCTCGGCATCCTTGCGGTTCTTCTCTGCCATCTTCTCATACTGGTCACGCATCTCGTTGAGGATGCGGCTCAGGTCCACGCCTGGGGCAGCGTCCATCTCCACATTGATCTCACCACCCACCTGGCCTCGCAGGGCGTTCATCTCCTCCTCGTGGTTCTTCTTCAGGTAG",
        "LCN2": "ATGGTGTTCGGGCTGGTGCGGCAGCTGGCGGCACCTGTGCACTCAGCCGTCGATACACTGGTCGATTGGGACAGGGAAGACGATGTGGTTTTCAGGGAGGCCCAGAGATTTGGAGAAGCGGATGAAGTTCTCCTTTAG",
        "MTAP": "ATGGATGAGCAACTGTGCCCGGCCAACAAAAAGCATTTCTATAAAATAACAATAAAAAAAGAGTTAAAAAAGCAAGCTTCCACTTTATTAAAAAGCAAAAATTACCAACTCGTATAA"
    }
    input_sequences = {}
    input_sequences = readingFile()
    encoded_seq, encoded_orf = dataEncoding(input=input_sequences, orf=Orf)
    counts(input_sequences, start, end)


main()
