from src.parser import NeedlemanWunschParser
from src.sequence_alignment import SequenceAlignment

seq1, seq2, n, input_filepath, gap_penalty, output_filename = NeedlemanWunschParser().get_parameters()

sequence_alignment = SequenceAlignment(seq1, seq2, input_filepath, strategy='global', gap_penalty=gap_penalty)
sequence_alignment.find_optimal_alignments(n, output_filename)
