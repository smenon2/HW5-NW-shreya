# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("../data/test_seq1.fa")
    seq2, _ = read_fasta("../data/test_seq2.fa")

    g = NeedlemanWunsch(sub_matrix_file='../substitution_matrices/BLOSUM62.mat', gap_open=-10, gap_extend=-1)
    g.align(seq1, seq2)


    test_align_matrix = np.array([[  0., -np.inf, -np.inf, -np.inf],
       [-np.inf,   5., -11., -13.],
       [-np.inf, -12.,   4.,  -8.],
       [-np.inf, -12.,  -1.,   5.],
       [-np.inf, -14.,  -6.,   4.]])

    test_gapA = np.array([[-10., -np.inf, -np.inf, -np.inf],
       [-11., -12.,  -6.,  -7.],
       [-12., -13., -14.,  -7.],
       [-13., -14., -15., -12.],
       [-14., -15., -16., -17.]])

    test_gapB = np.array([[-10., -11., -12., -13.],
       [-np.inf, -12., -13., -14.],
       [-np.inf,  -6., -14., -15.],
       [-np.inf,  -7.,  -7., -16.],
       [-np.inf,  -8.,  -8.,  -6.]])

    assert np.array_equal(g._align_matrix,test_align_matrix), "Align matrices don't match"
    assert np.array_equal(g._gapA_matrix, test_gapA), "Gap A matrices don't match"
    assert np.array_equal(g._gapB_matrix, test_gapB), "Gap B matrices don't match "
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("../data/test_seq3.fa")
    seq4, _ = read_fasta("../data/test_seq4.fa")

    g = NeedlemanWunsch(sub_matrix_file='../substitution_matrices/BLOSUM62.mat', gap_open=-10, gap_extend=-1)
    score, align3, align4 = g.align(seq3, seq4)

    assert score == 17, "score is wrong"
    assert align3 == 'MAVHQLIRRP', 'seq3 alignment wrong'
    assert align4 == 'M___QLIRHP', 'seq4 alignment wrong'




