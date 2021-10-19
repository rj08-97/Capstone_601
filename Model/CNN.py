from Bio import SeqIO, pairwise2
import pandas as pd
import numpy as np
from Bio.pairwise2 import format_alignment
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, MinMaxScaler
from sklearn.cluster import KMeans


def readingFile():
    seq = SeqIO.parse(open("/Users/rakshandajha/Documents/Capstone_601/DNA_Seq/GRCh38_latest_rna.fasta"), 'fasta')
    fasta_name = []
    fasta_sequence = []
    for s in seq:
        fasta_name.append(s.id)
        fasta_sequence.append(str(s.seq))
    return fasta_name, fasta_sequence


def fit(encoded, label):
    preprocessor = Pipeline(
        [
            ("scaler", MinMaxScaler()),
        ]
    )
    clusterer = Pipeline(
        [
            (
                "kmeans",
                KMeans(
                    n_clusters=3,
                    init="k-means++",
                    n_init=50,
                    max_iter=20,
                    random_state=42,
                ),
            ),
        ]
    )
    X = encoded.to_numpy()
    pipe = Pipeline(
        [
            ("preprocessor", preprocessor),
            ("clusterer", clusterer)
        ]
    )
    pipe.fit(X, label).transform()
    pipe.predict(X)


def dataEncoding(input):
    sequence = []
    mapping = dict(zip("ACGT", range(4)))
    for index in range(len(input)):
        sequence = [mapping[i] for i in input[index]]

    sequence = np.eye(4)[sequence]
    return sequence


def pairwise(sequence):
    sequence = list(sequence.values())
    align = []
    score = []
    l_score = []
    local_align = []
    for i in range(0, len(sequence) - 1):
        align = pairwise2.align.globalxx(sequence[i], sequence[i + 1])
        seq_length = min(len(sequence[i]), len(sequence[i + 1]))
        matches = align[0][2]
        score = (matches / seq_length) * 100
        print("The sequence similarity score for global alignment ", score)
        local_align = pairwise2.align.localxx(sequence[i], sequence[i + 1])
        local_matches = local_align[0][2]
        l_score = (local_matches / seq_length) * 100
        print("The sequence similarity score for local alignment ", l_score)
    print("Needleman-Wunsch Global Algorithm")
    for index in range(len(align)):
        print(format_alignment(*align[index]))
    print("Smith-Waterman Local Algorithm")
    for index in range(len(local_align)):
        print(format_alignment(*local_align[index]))
    return np.array(align), np.array(local_align), score, l_score


def validation(sequence):
    dna_char = "ACGT"
    for i in range(len(sequence)):
        for letter in sequence[i]:
            if letter not in dna_char:
                print("Not a valid sequence")
                return False
        return True


def preprocess(data):
    df = data.iloc[:, 1:]
    label = data.iloc[:, :1]
    return df, label


def main():
    mRNA = pd.read_excel("/Users/rakshandajha/Documents/Capstone_601/DNA_Seq/overlap.xlsx")
    name, sequence = readingFile()
    seq = sequence[0:100]
    validation(sequence)
    df, label = preprocess(mRNA)
    encoded = dataEncoding(seq)
    print(encoded)
    fit(df, label)


main()
