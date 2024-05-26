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
            identifier, sequence, quality = '', '', ''
            for line in file:
                line = line.strip()
                if file_type == 'fasta':
                    if line.startswith('>'):
                        if identifier:
                            sequences[identifier] = sequence
                            print(f"DEBUG: Added sequence {identifier} with length {len(sequence)}")
                        identifier = line[1:]
                        sequence = ''
                    else:
                        sequence += line
                elif file_type == 'fastq':
                    if line.startswith('@') and not identifier:
                        identifier = line[1:]
                    elif line.startswith('+') and identifier:
                        continue
                    elif not sequence and not quality:
                        sequence += line
                    elif sequence and not quality:
                        quality += line
                    if quality and len(quality) >= len(sequence):
                        sequences[identifier] = (sequence, quality)
                        print(f"DEBUG: Added sequence {identifier} with length {len(sequence)} and quality length {len(quality)}")
                        identifier, sequence, quality = '', '', ''
            if identifier:
                if file_type == 'fasta':
                    sequences[identifier] = sequence
                    print(f"DEBUG: Added final sequence {identifier} with length {len(sequence)}")
                else:
                    sequences[identifier] = (sequence, quality)
                    print(f"DEBUG: Added final sequence {identifier} with length {len(seqeunce)} and quality length {len(quality)}")
        return sequences
    except Exception as e:
        sys.stderr.write(f"ERROR: Failed to parse {file_type.upper()} file {file_path}: {str(e)}\n")
        sys.exit(1)
        
def align_references(queries, reference_id, reference_seq):
    aligner = align_it(reference_seq)
    results = []
    for query_id, query_info in queries.items():
        query_seq, query_qual = query_info
        match_positions = aligner.search(query_seq)
        if match_positions:
            for pos in match_positions:
                results.append(f"{query_id}\t0\t{reference_id}\t{pos+1}\t255\t{len(query_seq)}M\t*\t0\t0\t{query_seq}\t{query_qual}")
        else:
            results.append(f"{query_id}\t4\t*\t0\t0\t*\t*\t0\t0\t{query_seq}\t{query_qual}")
    return results

def main():
    parser = argparse.ArgumentParser(description="Align FASTQ reads against a reference genome.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ file containing reads.")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome in FASTA format.")
    
    args = parser.parse_args()
    references = parse_sequences(args.reference, 'fasta')
    queries = parse_sequences(args.input, 'fastq')  # This now returns a dict with (sequence, quality)

    print("@HD\tVN:1.6\tSO:unsorted")
    for ref_id in references:
        print(f"@SQ\tSN:{ref_id}\tLN:{len(references[ref_id])}")

    with ThreadPoolExecutor() as executor:
        futures = {executor.submit(align_references, queries, ref_id, references[ref_id]): ref_id for ref_id in references}
        for future in as_completed(futures):
            for result in future.result():
                print(result)

if __name__ == "__main__":
    main()

