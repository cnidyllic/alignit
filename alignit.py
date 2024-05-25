import sys

class alignit:
    def __init__(self, reference_seq, k=20):
        self.reference_seq = reference_seq
        self.k = k
        self.index = self.build_index()

    def build_index(self):
        """Builds a hash-based index for k-mers from the reference sequence with overlap handling."""
        index = {}
        step = max(1, self.k // 4)  # Reduce the step to handle overlapping k-mers efficiently
        for i in range(0, len(self.reference_seq) - self.k + 1, step):
            kmer = self.reference_seq[i:i+self.k]
            index.setdefault(kmer, []).append(i)
        return index

    def search(self, query):
        """Searches for the query in the reference sequence using the hash-based index."""
        matches = set()
        for i in range(len(query) - self.k + 1):
            kmer = query[i:i+self.k]
            if kmer in self.index:
                for pos in self.index[kmer]:
                    end_pos = pos + len(query)
                    if end_pos <= len(self.reference_seq) and query == self.reference_seq[pos:end_pos]:
                        matches.add(pos)
        return sorted(list(matches))

def parse_fasta(fasta_file):
    """Parses a FASTA file and returns a dictionary of sequences."""
    try:
        sequences = {}
        with open(fasta_file, 'r') as file:
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
        handle_error(f"Failed to parse FASTA file {fasta_file}: {str(e)}")

def parse_fastq(fastq_file):
    """Parses a FASTQ file and returns a list of (identifier, sequence) tuples."""
    try:
        sequences = []
        with open(fastq_file, 'r') as file:
            while True:
                identifier = file.readline().strip()
                if not identifier:
                    break
                sequence = file.readline().strip()
                file.readline()  # Plus line (ignored)
                file.readline()  # Quality line (ignored)
                sequences.append((identifier[1:], sequence))  # Remove '@' character
        return sequences
    except Exception as e:
        handle_error(f"Failed to parse FASTQ file {fastq_file}: {str(e)}")

def handle_error(message):
    """Handle errors by logging and exiting, with improved logging for diagnostics."""
    sys.stderr.write(f"ERROR: {message}\n")
    sys.exit(1)

def main():
    try:
        reference_path = "example-files/test_reference.fa"
        query_path = "example-files/test_queries.fastq"

        references = parse_fasta(reference_path)
        queries = parse_fastq(query_path)

        for ref_id, ref_seq in references.items():
            print(f"Aligning to reference: {ref_id}")
            aligner = alignit(ref_seq)
            for query_id, query_seq in queries:
                match_positions = aligner.search(query_seq)
                if match_positions:
                    print(f"{query_id} found in {ref_id} at positions: {match_positions}")
                else:
                    print(f"{query_id} not found in {ref_id}.")
    except Exception as e:
        handle_error(f"An unexpected error occurred: {str(e)}")

if __name__ == "__main__":
    main()
