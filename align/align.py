# Importing Dependencies
import numpy as np
from typing import Tuple


# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """

    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file)  # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB

        # TODO: Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing

        # create matrices and initialize to negative infinity
        # align_matrix keeps track of scores if aligned sequences
        # gapA keeps track of score if gap in seqA
        # gapB keeps track of score if gap in seqB
        # back keeps track for backtracing

        # You need to initialize to infinity because there are values should be kept empty like first column in gapB
        # (zero/one initialize won't work)

        self._align_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapA_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapB_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._back = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf

        # TODO: Implement global alignment here
        # Set the starting node (0,0) to zero in the align matrices
        self._align_matrix[0][0] = 0

        # Initialize first row and and column of gapA and gapB as if all gaps
        for i in range(0, len(seqA) + 1):
            self._gapA_matrix[i][0] = self.gap_open + (self.gap_extend * i)

        for j in range(0, len(seqB) + 1):
            self._gapB_matrix[0][j] = self.gap_open + (self.gap_extend * j)

        # Now loop through sequence A and sequence B to fill in all the matrices
        for i in range(1, len(seqA) + 1):
            for j in range(1, len(seqB) + 1):

                # Get base at position
                A = self._seqA[i - 1]
                B = self._seqB[j - 1]

                # Get the score from the scoring matrices of aligning A and B
                p_score = self.sub_dict[(A, B)]

                # Now get the values of the three possibilities: match, gap in A, gap in B
                # Affine gap penalties

                # First calculate the score for matching:

                # Fill in alignment matrix
                gapinA = self._gapA_matrix[i - 1][j - 1]
                gapinB = self._gapB_matrix[i - 1][j - 1]
                p = self._align_matrix[i - 1][j - 1]
                self._align_matrix[i][j] = p_score + max(p, gapinA, gapinB)

                # Fill in gapA_matrix: extend the gap
                gapA1 = self.gap_extend + self._gapA_matrix[i][j - 1]
                gapA2 = self.gap_open + self.gap_extend + self._gapB_matrix[i][j - 1]
                gapA3 = self.gap_open + self.gap_extend + self._align_matrix[i][j - 1]
                self._gapA_matrix[i][j] = max(gapA1, gapA2, gapA3)

                # Fill in gapB_matrix
                gapB1 = self.gap_extend + self._gapB_matrix[i - 1][j]
                gapB2 = self.gap_open + self.gap_extend + self._gapA_matrix[i - 1][j]
                gapB3 = self.gap_open + self.gap_extend + self._align_matrix[i - 1][j]
                self._gapB_matrix[i][j] = max(gapB1, gapB2, gapB3)

                # Get the best score from the three possibilities
                max_score = max(self._align_matrix[i][j], self._gapA_matrix[i][j], self._gapB_matrix[i][j])

                # We are going to store the pointers in the back matrix
                # Match is diag, 0
                # gap in A is up, 1
                # gap in B is down, 2

                if max_score == self._align_matrix[i][j]:
                    self._back[i][j] = 0
                elif max_score == self._gapA_matrix[i][j]:
                    self._back[i][j] = -1
                else:
                    self._back[i][j] = 1

        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Get score - which is stored in align matrix
        m = len(self._seqA)
        n = len(self._seqB)
        self.alignment_score = max(self._align_matrix[m][n], self._gapB_matrix[m][n], self._gapA_matrix[m][n])

        # Trying to loop backwards thru seqA and seqB:
        while m > 0 or n > 0:
            pointer = self._back[m][n]
            # If pointer is 0, that means that the best score was a match (or mismatch, just not gap)
            if pointer == 0:
                self.seqA_align = self._seqA[m - 1] + self.seqA_align
                self.seqB_align = self._seqB[n - 1] + self.seqB_align
                m = m - 1
                n = n - 1

            # If pointer is -1, then there was a gap in A
            elif pointer == -1:
                self.seqA_align = "_" + self.seqA_align
                self.seqB_align = self._seqB[n - 1] + self.seqB_align
                n = n - 1

            # If pointer is 1, then there was a gap in B
            else:
                self.seqA_align = self._seqA[m - 1] + self.seqA_align
                self.seqB_align = "_" + self.seqB_align
                m = m - 1

        return self.alignment_score, self.seqA_align, self.seqB_align


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
