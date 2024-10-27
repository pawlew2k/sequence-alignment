# Needleman–Wunsch algorithm: https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
# Smith–Waterman algorithm: https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm


import os
from collections import defaultdict
from dataclasses import dataclass

import numpy as np
import pandas as pd
from typing import List, Literal, Tuple

from src.directional_cell import DirectionalCell


@dataclass
class SequenceAlignment:
    """A class implementing the Needleman–Wunsch and Smith–Waterman algorithms to perform alignment of two DNA sequences

    Needleman–Wunsch == strategy='global'
    Smith–Waterman == strategy='local'

    Args:
        seq1: The first nucleotide sequence
        seq2: The second nucleotide sequence
        input_filepath: File path to the substitution matrix in CSV format (see template substitution_matrix.csv)
        strategy: Global or local alignment strategy
        gap_penalty: The value of the constant penalty for inserting a gap ('-')
    """

    def __init__(self, seq1: str, seq2: str, input_filepath: str,
                 strategy: Literal['global', 'local'], gap_penalty: int = -2):
        self.seq1 = seq1
        self.seq2 = seq2
        self.strategy = strategy
        self.gap_penalty = gap_penalty

        self.alignments = self._load_substitution_matrix(input_filepath)

    def find_optimal_alignments(self, n: int, output_filename: str):
        """Returns n optimal alignments based on the best score"""
        self._create_score_and_directional_matrices()
        solution_sequences, score = self._traceback(n)
        self._print_and_save_solutions(solution_sequences, score, n, output_filename)

    @staticmethod
    def _load_substitution_matrix(file_path: str):
        # Load the substitution matrix from a CSV file
        full_matrix = pd.read_csv(file_path, header=None)
        seq1 = ''.join(full_matrix.iloc[1:, 0])
        seq2 = ''.join(full_matrix.iloc[0, 1:])
        substitution_matrix = np.asarray(full_matrix.iloc[1:, 1:], dtype=int)

        # Score for all possible pairwise alignments of each nucleotide with every other nucleotide
        alignments = defaultdict(lambda: defaultdict(lambda x: 0))
        for i, nucleotide1 in enumerate(seq1):
            for j, nucleotide2 in enumerate(seq2):
                alignments[nucleotide1][nucleotide2] = substitution_matrix[i][j]
        return alignments

    def _create_score_and_directional_matrices(self):
        seq1, seq2 = self.seq1, self.seq2

        self.score_matrix = np.zeros(shape=(len(seq1) + 1, len(seq2) + 1), dtype=int)
        self.directional_matrix = np.full(shape=(len(seq1) + 1, len(seq2) + 1), fill_value=DirectionalCell())

        # Filling the edges with the score and directional matrices
        for i in range(1, len(seq1) + 1):
            self.score_matrix[i, 0] = i * self.gap_penalty
            self.directional_matrix[i, 0] = DirectionalCell(up=True)
        for i in range(1, len(seq2) + 1):
            self.score_matrix[0, i] = i * self.gap_penalty
            self.directional_matrix[0, i] = DirectionalCell(left=True)

        # Dynamic programming matrix filling
        for i, nucleotide1 in enumerate(seq1):
            for j, nucleotide2 in enumerate(seq2):
                # Scores from possible three directions
                diagonal_score = self.score_matrix[i, j] + self.alignments[nucleotide1][nucleotide2]
                up_score = self.score_matrix[i, j + 1] + self.gap_penalty
                left_score = self.score_matrix[i + 1, j] + self.gap_penalty

                # Save the best score from the directions
                best_score = np.max([up_score, diagonal_score, left_score])
                if self.strategy == 'local':  # Minimum 0 in local alignment
                    best_score = max(best_score, 0)
                self.score_matrix[i + 1, j + 1] = best_score

                # Create an object representing possible paths for backtracking
                directional_cell = DirectionalCell(up=(up_score == best_score),
                                                   diagonal=(diagonal_score == best_score),
                                                   left=(left_score == best_score))
                self.directional_matrix[i + 1, j + 1] = directional_cell

    def _traceback(self, n) -> Tuple[list, int]:
        # Create solution states representing the current cell with trace [posY, posX, seq1, seq2]
        match self.strategy:
            case 'global':  # From the bottom right cell
                max_score = self.score_matrix[-1, -1]
                solution_states = [list(self.directional_matrix[1:, 1:].shape) + ['', '']]
            case 'local':  # All cells with the maximum score
                max_score = np.max(self.score_matrix)
                solution_states = [cell + ['', ''] for cell in np.argwhere(self.score_matrix == max_score).tolist()]
            case _:
                raise AttributeError("Strategy must be 'global' or 'local'.")

        # Go deeper into the path you last traveled and add new possible solutions
        solution_sequences = set()
        while solution_states and len(solution_sequences) < n:
            posY, posX, seq1, seq2 = solution_states.pop()
            cell = self.directional_matrix[posY][posX]

            # Finish travel if no more possible directions or value of cell equals to 0 in local alignment
            if cell.num_directions == 0 or (self.strategy == 'local' and self.score_matrix[posY, posX] == 0):
                solution_sequences.add((seq1, seq2))

            # Add new possible solutions in three directions if there are any
            # Update cell position and sequences adding a single nucleotide or gap '-'
            else:
                up, diagonal, left = cell.get_directions()
                if left:
                    solution_states.append([posY, posX - 1, seq1 + '-', seq2 + self.seq2[-posX]])
                if up:
                    solution_states.append([posY - 1, posX, seq1 + self.seq1[-posY], seq2 + '-'])
                if diagonal:
                    solution_states.append([posY - 1, posX - 1, seq1 + self.seq1[-posY], seq2 + self.seq2[-posX]])

        return sorted(solution_sequences, reverse=True), max_score

    def _print_and_save_solutions(self, solution_sequences: List[Tuple[str, str]], score: int, n: int, output_filename: str, output_dir='output'):
        # Print the solutions
        text = ''
        for i in range(min(n, len(solution_sequences))):
            seq1, seq2 = solution_sequences[i]
            text += (f'{self.strategy.capitalize()} alignment no. {i + 1}:\n'
                     f'{seq1}\n'
                     f'{seq2}\n'
                     f'Score: {score}\n\n')
        text = text[:-2]
        print(text)

        # Create directory if needed
        if output_dir not in os.listdir():
            os.mkdir(output_dir)

        # Save the solution
        with open(os.path.join(output_dir, output_filename), 'w') as f:
            f.write(text)
