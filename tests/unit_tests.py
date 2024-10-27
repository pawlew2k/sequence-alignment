import argparse
import unittest
from unittest.mock import patch

import pandas as pd

from src.directional_cell import DirectionalCell
from src.sequence_alignment import SequenceAlignment
from src.parser import NeedlemanWunschParser, SmitchWatermanParser


class TestDirectionalCell(unittest.TestCase):
    def test_direction_initialization(self):
        # Test with all directions True
        cell = DirectionalCell(up=True, diagonal=True, left=True)
        self.assertEqual(cell.num_directions, 3)
        self.assertEqual(cell.get_directions(), (True, True, True))

        # Test with only 'up' direction True
        cell = DirectionalCell(up=True)
        self.assertEqual(cell.num_directions, 1)
        self.assertEqual(cell.get_directions(), (True, False, False))

        # Test with no directions True
        cell = DirectionalCell()
        self.assertEqual(cell.num_directions, 0)
        self.assertEqual(cell.get_directions(), (False, False, False))


class TestSequenceAlignment(unittest.TestCase):
    # Mock CSV data for substitution matrix
    mock_data = [["", "A", "C", "G", "T"],
                 ["A", 1, -1, -2, -1],
                 ["C", -1, 1, -1, -2],
                 ["G", -2, -1, 1, -1],
                 ["T", -1, -2, -1, 1]]

    @patch('sequence_alignment.pd.read_csv')
    def test_load_substitution_matrix(self, mock_read_csv):
        mock_read_csv.return_value = pd.DataFrame(self.mock_data)
        sequence_alignment = SequenceAlignment(seq1="ACGT", seq2="TGCA", input_filepath="dummy.csv", strategy="global", gap_penalty=-2)
        alignments = sequence_alignment._load_substitution_matrix("dummy.csv")

        self.assertEqual(alignments['A']['A'], 1)
        self.assertEqual(alignments['A']['T'], -1)
        self.assertEqual(alignments['G']['C'], -1)

    @patch('sequence_alignment.pd.read_csv')
    def test_create_score_and_directional_matrices(self, mock_read_csv):
        mock_read_csv.return_value = pd.DataFrame(self.mock_data)

        alignment = SequenceAlignment(seq1="A", seq2="T", input_filepath="dummy.csv", strategy="global", gap_penalty=-2)
        alignment.alignments = {'A': {'A': 1, 'T': -1}, 'T': {'A': -1, 'T': 1}}
        alignment._create_score_and_directional_matrices()

        # Check score matrix shape and values
        self.assertEqual(alignment.score_matrix.shape, (2, 2))
        self.assertEqual(alignment.score_matrix[1, 1], -1)

        # Check directional matrix types
        self.assertIsInstance(alignment.directional_matrix[1, 1], DirectionalCell)

    @patch('sequence_alignment.pd.read_csv')
    def test_find_optimal_alignments_global(self, mock_read_csv):
        mock_read_csv.return_value = pd.DataFrame(self.mock_data)

        alignment = SequenceAlignment(seq1="A", seq2="A", input_filepath="dummy.csv", strategy="global", gap_penalty=-2)
        alignment.alignments = {'A': {'A': 1}}

        # Run optimal alignment search
        alignment._create_score_and_directional_matrices()
        solutions, score = alignment._traceback(1)

        self.assertEqual(len(solutions), 1)
        self.assertEqual(score, 1)
        self.assertEqual(solutions[0], ("A", "A"))


class TestParsers(unittest.TestCase):
    @patch('argparse.ArgumentParser.parse_args')
    def test_needleman_wunsch_parser(self, mock_args):
        mock_args.return_value = argparse.Namespace(seq1="ACGT", seq2="TGCA", n=1, substitution_matrix="matrix.csv",
                                                    gap_penalty=-2, output="output.txt")

        parser = NeedlemanWunschParser()
        params = parser.get_parameters()

        self.assertEqual(params, ("ACGT", "TGCA", 1, "matrix.csv", -2, "output.txt"))

    @patch('argparse.ArgumentParser.parse_args')
    def test_smith_waterman_parser(self, mock_args):
        mock_args.return_value = argparse.Namespace(seq1="ACGT", seq2="TGCA", n=2, substitution_matrix="matrix.csv",
                                                    gap_penalty=-1, output="output_sw.txt")

        parser = SmitchWatermanParser()
        params = parser.get_parameters()

        self.assertEqual(params, ("ACGT", "TGCA", 2, "matrix.csv", -1, "output_sw.txt"))


if __name__ == '__main__':
    unittest.main()
