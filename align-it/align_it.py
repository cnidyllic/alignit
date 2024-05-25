import sys
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

class align_it:
    def __init__(self, reference_seq, k=20):
        self.reference_seq = reference_seq
        self.k = self.determine_kmer_size(reference_seq)
        self.index = self.build_index()
        
    def determine_kmer_size(self, sequence):
        """ Dynamically adjust k-mer size based on GC content to optimize indexing. """
        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
        if gc_content < 0.4:
            return max(20, len(sequence) // 100)  # Increase k-mer size for high AT-content
        return 20
        
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
        return (kmer.count('G') + kmer.count('C')) / len(kmer) > 0.2  # Threshold for significance, can be adjusted

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
        
def parse_sequences(file_path, file_type):
    try:
        sequences = {}
        with open(file_path, 'r') as file:
            identifier, sequence = '', ''
            for line in file:
                line = line.strip()
                if line.startswith('>' if file_type == 'fasta' else '@'):
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
        sys.stderr.write(f"ERROR: Failed to parse {file_type.upper()} file {file_path}: {str(e)}\n")
        sys.exit(1)

def align_reads(ref_seq, queries):
    aligner = align_it(ref_seq)
    results = {}
    for query_id, query_seq in queries.items():
        match_positions = aligner.search(query_seq)
        if match_positions:
            results[query_id] = (True, match_positions)
        else:
            results[query_id] = (False, [])
    return results

def main():
    parser = argparse.ArgumentParser(description="Align FASTQ reads against a reference genome.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ file containing reads.")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome in FASTA format.")
    args = parser.parse_args()

    reference_sequences = parse_sequences(args.reference, 'fasta')
    queries = parse_sequences(args.input, 'fastq')
    total_reads = len(queries)

    results = {}
    with ThreadPoolExecutor() as executor:
        future_to_ref = {executor.submit(align_reads, ref_seq, queries): ref_id for ref_id, ref_seq in reference_sequences.items()}
        for future in as_completed(future_to_ref):
            ref_id = future_to_ref[future]
            results[ref_id] = future.result()

    reads_aligned = sum(1 for res in results.values() for status, positions in res.items() if status)
    if total_reads > 0:
        alignment_percentage = (reads_aligned / total_reads) * 100
        print(f"Percentage of reads aligned: {alignment_percentage:.2f}%")
    else:
        print("No reads to align.")

if __name__ == "__main__":
    main()
