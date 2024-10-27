from src.parser import SmitchWatermanParser
from src.sequence_alignment import SequenceAlignment

seq1, seq2, n, input_filepath, gap_penalty, output_filename = SmitchWatermanParser().get_parameters()

sequence_alignment = SequenceAlignment(seq1, seq2, input_filepath, strategy='local', gap_penalty=gap_penalty)
sequence_alignment.find_optimal_alignments(n, output_filename)
