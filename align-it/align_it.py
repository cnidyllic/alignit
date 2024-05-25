import sys
import argparse

class align_it:
    def __init__(self, reference_seq, k=20):
        self.reference_seq = reference_seq
        self.k = k
        self.index = self.build_index()

    def build_index(self):
        """ Builds a hash-based index for k-mers from the reference sequence with optimized handling for AT-rich sequences. """
        index = {}
        step = max(1, self.k // 4)  # Adjusted step for overlap, could be increased based on AT-content analysis
        for i in range(0, len(self.reference_seq) - self.k + 1, step):
            kmer = self.reference_seq[i:i+self.k]
            if self.is_significant_kmer(kmer):  # Check if the k-mer is significant enough to be indexed
                index.setdefault(kmer, []).append(i)
        return index

    def is_significant_kmer(self, kmer):
        """ Determine if a k-mer is significant based on its GC content in high AT-context. """
        gc_content = (kmer.count('G') + kmer.count('C')) / len(kmer)
        return gc_content > 0.2  # Threshold for significance, can be adjusted

    def search(self, query):
        """ Searches for the query in the reference sequence using the hash-based index. """
        matches = set()
        for i in range(len(query) - self.k + 1):
            kmer = query[i:i+self.k]
            if kmer in self.index:
                for pos in self.index[kmer]:
                    end_pos = pos + len(query)
                    if end_pos <= len(self.reference_seq) and query == self.reference_seq[pos:end_pos]:
                        matches.add(pos)
        return sorted(list(matches))

def parse_fasta(file_path):
    """ Parses a FASTA file and returns a dictionary of sequences. """
    try:
        sequences = {}
        with open(file_path, 'r') as file:
            identifier, sequence = '', ''
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if identifier:
                        sequences[identifier] = sequence
                    identifier = line[1:]
                    sequence = ''
                else:
                    sequence += line
            if identifier:
                sequences[identifier] = sequence
        return sequences
    except Exception as e:
        sys.stderr.write(f"ERROR: Failed to parse FASTA file {file_path}: {str(e)}\n")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Align FASTQ reads against a reference genome.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ file containing reads.")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome in FASTA format.")
    
    args = parser.parse_args()

    reference_sequences = parse_fasta(args.reference)
    queries = parse_fastq(args.input)

    for ref_id, ref_seq in references.items():
        print(f"Aligning to reference: {ref_id}")
        aligner = align_it(ref_seq)
        for query_id, query_seq in queries.items():
            match_positions = aligner.search(query_seq)
            if match_positions:
                print(f"{query_id} found in {ref_id} at positions: {match_positions}")
            else:
                print(f"{query_id} not found in {ref_id}.")

if __name__ == "__main__":
    main()
