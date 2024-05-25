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

def parse_fastq(fastq_file):
    """
    Parses a FASTQ file and returns a list of tuples containing the read identifier,
    sequence, and quality string.
    
    Parameters:
        fastq_file (str): Path to the FASTQ file.
    
    Returns:
        list of tuples: List containing tuples of (identifier, sequence, quality).
    """
    queries = []
    try:
        with open(fastq_file, 'r') as file:
            while True:
                identifier = file.readline().strip()
                if not identifier:
                    break  # End of file
                if not identifier.startswith('@'):
                    raise ValueError("Malformed FASTQ file: Expected '@' at the start of the identifier line.")
                
                sequence = file.readline().strip()
                separator = file.readline().strip()
                if not separator.startswith('+'):
                    raise ValueError("Malformed FASTQ file: Expected '+' at the start of the separator line.")
                
                quality = file.readline().strip()
                if len(sequence) != len(quality):
                    raise ValueError("Malformed FASTQ file: Sequence and quality scores length mismatch.")
                
                queries.append((identifier[1:], sequence, quality))  # Strip '@' from the identifier

    except IOError as e:
        sys.stderr.write(f"ERROR: Unable to open or read the FASTQ file {fastq_file}: {str(e)}\n")
        sys.exit(1)
    except ValueError as e:
        sys.stderr.write(f"ERROR: {str(e)}\n")
        sys.exit(1)

    return queries

def main():
    parser = argparse.ArgumentParser(description="Align FASTQ reads against a reference genome.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ file containing reads.")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome in FASTA format.")
    
    args = parser.parse_args()

    reference_sequences = parse_fasta(args.reference)
    queries = parse_fastq(args.input)

    total_reads = len(queries)
    reads_aligned = 0

    for ref_id, ref_seq in reference_sequences.items():
        print(f"Aligning to reference: {ref_id}")
        aligner = aligh_it(ref_seq)
        for query_id, query_seq in queries.items():
            match_positions = aligner.search(query_seq)
            if match_positions:
                print(f"{query_id} found in {ref_id} at positions: {match_positions}")
                reads_aligned += 1
            else:
                print(f"{query_id} not found in {ref_id}.")

    if total_reads > 0:
        alignment_percentage = (reads_aligned / total_reads) * 100
        print(f"Percentage of reads aligned: {alignment_percentage:.2f}%")
    else:
        print("No reads to align.")

if __name__ == "__main__":
    main()
