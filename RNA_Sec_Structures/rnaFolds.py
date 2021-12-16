"""Author - Rakshanda Jha"""
"""The following code produces the RNA secondary structures in dot-bracket notation along with the mapping of the 
binding sites. """

import os
from Bio import SeqIO
import rna_tools as rna
from rna_tools import Seq
import subprocess
from rna_tools.SecondaryStructure import draw_ss, parse_vienna_to_pairs
from rna_tools.Seq import RNASequence


def readingFile():
    """The function reads in the sequences in fasta format and returns it as a dictionary"""

    global sequences
    sequence = SeqIO.parse(open("../DNA_Seq/AVIL_mRNA.fasta"), 'fasta')
    for fasta in sequence:
        name, sequences = fasta.id, str(fasta.seq)
        # new_sequence = some_function(sequence)


    return sequences


def rnaStructure(input):
    """The function takes in the sequence as an input and returns the RNA secondary structure in a dot-bracket
    notation along with the binding sites positions. """
    ss = []
    motif = ['CAGUGA', 'GCGC', 'ACUCU', 'UGUU', 'UUGUG', 'UUGUU', 'AUGUG', 'CAUCG', 'GUGUG']
    for seq in input.values():
        seq = RNASequence(seq)
        # print output
        print(seq.predict_ss('centroid_fold'))
        print("RNA fold")
        s = seq.predict_ss(method="RNAfold")
        # print(d)
        ss = s.split("\n")
        print("RNA Structure" + "\n" + s)
    position = []
    structure = {}
    i = 0
    # mapping
    for seq in ss:
        for mot in motif:
            if mot in seq:
                position.append(seq.index(mot))
                structure = {mot, seq[i]}
                # print(structure)
                print(position)
        i = i + 1
    return structure


def main():
    input_sequences = readingFile()
    rnaStructure(input_sequences)


main()
