import argparse
from typing import Tuple


class Parser:
    def __init__(self):
        self.parser = None

    def _add_arguments(self):
        self.parser.add_argument('--seq1', type=str, required=True, help='first sequence of nucleotides')
        self.parser.add_argument('--seq2', type=str, required=True, help='second sequence of nucleotides')
        self.parser.add_argument('-n', type=int, required=True, help='count of maximum optimal alignments')
        self.parser.add_argument('-s', '--substitution_matrix', type=str, required=True,
                                 help='filepath to the substitution matrix in CSV format')
        self.parser.add_argument('-g', '--gap_penalty', type=int, default=-2,
                                 help="value of constant penalty for insertion a gap ('-') incurs")
        self.parser.add_argument('-o', '--output', type=str, required=True,
                                 help='output filename, file will be saved in \'/output\' directory')

    def get_parameters(self) -> tuple:
        return tuple(vars(self.parser.parse_args()).values())


class NeedlemanWunschParser(Parser):
    def __init__(self):
        super().__init__()
        self.parser = argparse.ArgumentParser(description='Needleman-Wunsch (NW) algorithm for global alignment')
        self._add_arguments()


class SmitchWatermanParser(Parser):
    def __init__(self):
        super().__init__()
        self.parser = argparse.ArgumentParser(description='Smith-Waterman (SW) algorithm for local alignment')
        self._add_arguments()
