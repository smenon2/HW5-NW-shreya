# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch
import numpy as np


def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    g = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', gap_open=-10, gap_extend=-1)
    ggScore, gg1, gg2 = g.align(gg_seq, hs_seq)
    mmScore, mm1, mm2 = g.align(mm_seq, hs_seq)
    brScore, br1, br2 = g.align(br_seq, hs_seq)
    ttScore, tt1, tt2 = g.align(tt_seq, hs_seq)

    species = np.array(['Gallus gallus', 'Mus musculus', 'Balaeniceps rex', 'tursiops truncatus'])
    scores = np.array([ggScore, mmScore, brScore, ttScore])
    idx = np.argsort(scores)[::-1]
    print(species[idx])

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    for i in range(0, len(scores)):
        print(species[idx][i], scores[idx][i])


if __name__ == "__main__":
    main()
