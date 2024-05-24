from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import sys

class SequenceAligner:
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file

    def load_sequences(self):
        """Load sequences from a FASTA file."""
        try:
            sequences = list(SeqIO.parse(self.input_file, "fasta"))
            print(f"Loaded {len(sequences)} sequences from {self.input_file}")
            return sequences
        except Exception as e:
            print(f"An error occurred while loading sequences: {e}")
            sys.exit(1)

    def align_sequences(self, sequences):
        """Align sequences using a simple pairwise method."""
        print("Aligning sequences with SequenceAligner Tool...")

        # Check if there are at least two sequences to align
        if len(sequences) < 2:
            raise ValueError("At least two sequences are needed for alignment")

        # Aligning using pairwise2 
        alignments = pairwise2.align.globalxx(sequences[0].seq, sequences[1].seq)

        # Select the best alignment (you could also choose to handle multiple alignments)
        best_alignment = alignments[0]
        aligned_seq1 = format_alignment(*best_alignment).split('\n')[0]
        aligned_seq2 = format_alignment(*best_alignment).split('\n')[1]

        # Create a MultipleSeqAlignment object to return
        aligned_sequences = MultipleSeqAlignment([
            SeqIO.SeqRecord(Seq(aligned_seq1.replace(' ', '-')), id=sequences[0].id),
            SeqIO.SeqRecord(Seq(aligned_seq2.replace(' ', '-')), id=sequences[1].id)
        ])

        return aligned_sequences

    def save_aligned_sequences(self, aligned_sequences):
        """Save the aligned sequences to an output file."""
        try:
            AlignIO.write(aligned_sequences, self.output_file, "fasta")
            print(f"Aligned sequences saved to {self.output_file}")
        except Exception as e:
            print(f"An error occurred while saving the aligned sequences: {e}")
            sys.exit(1)

    def run_alignment(self):
        """Run the full alignment process."""
        sequences = self.load_sequences()
        aligned_sequences = self.align_sequences(sequences)
        self.save_aligned_sequences(aligned_sequences)

# Usage 
aligner = SequenceAligner("path_to_input.fasta", "path_to_output.fasta")
aligner.run_alignment()
