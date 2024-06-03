import sys
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from math import log2, ceil

class align_it:
    def __init__(self, reference_seq, k_override=None, significance_threshold=0.2, entropy_threshold=1.5):
        self.reference_seq = reference_seq
        self.significance_threshold = significance_threshold
        self.entropy_threshold = entropy_threshold
        self.k = koverride or 20
        self.index = self.build_index(self.reference_seq, self.k)

    def calculate_entropy(self, kmer):
        """ Calculates the Shannon entropy of the k-mer. """
        freq = {x: kmer.count(x) / len(kmer) for x in set(kmer)}
        entropy = -sum(p * log2(p) for p in freq.values() if p > 0)
        return entropy
        
    def adjust_kmer_size(self, sequence):
        """ Adjusts k-mer size based on GC content unless overridden by command-line. """
        if self.k is not None:
            return self.k  # Return overridden k-mer size if set
        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
        if gc_content < 0.4:
            return max(20, len(sequence) // 100)
        elif gc_content > 0.6:
            return max(15, len(sequence) // 150)
        return 20
        
    def build_index(self, sequence, k):
        """ Builds a hash-based index for k-mers from the reference sequence with optimized handling for AT-rich sequences. """
        index = {}
        step = max(1, k // 4)  # Adjusted step for overlap, could be increased based on AT-content analysis
        for i in range(0, len(sequence) - k + 1, step):
            kmer = sequence[i:i+k]
            if self.is_significant_kmer(kmer): 
                index.setdefault(kmer, []).append(i)
        return index

    def is_significant_kmer(self, kmer):
        """ Determine if a k-mer is significant based on its GC content and entropy. """
        gc_content_significant = (kmer.count('G') + kmer.count('C')) / len(kmer) > self.significance_threshold
        entropy_significant = self.calculate_entropy(kmer) > self.entropy_threshold
        return gc_content_significant and entropy_significant

    def search(self, query):
        """ Searches for the query in the reference sequence using a dynamic k-mer index. """
        k = self.adjust_kmer_size(query)
        self.index = self.build_index(self.reference_seq, k)
        matches = set()
        for i in range(len(query) - k + 1):
            kmer = query[i:i+k]
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
                        identifier, sequence, quality = '', '', ''
            if identifier:
                if file_type == 'fasta':
                    sequences[identifier] = sequence
                else:
                    sequences[identifier] = (sequence, quality)
        return sequences
    except Exception as e:
        sys.stderr.write(f"ERROR: Failed to parse {file_type.upper()} file {file_path}: {str(e)}\n")
        sys.exit(1)
        
def align_block(queries_block, references, k_override, significance_threshold, entropy_threshold):
    results = []
    for reference_id, reference_seq in references.items():
        aligner = align_it(reference_seq, k_override, significance_threshold, entropy_threshold)
        for query_id, query_info in queries_block.items():
            query_seq, query_qual = query_info
            match_positions = aligner.search(query_seq)
            if match_positions:
                for pos in match_positions:
                    results.append(f"{query_id}\t0\t{reference_id}\t{pos+1}\t255\t{len(query_seq)}M\t*\t0\t0\t{query_seq}\t{query_qual}")
            else:
                results.append(f"{query_id}\t4\t*\t0\t0\t*\t*\t0\t0\t{query_seq}\t{query_qual}")
    return results

def distribute_queries(queries, num_blocks):
    """ Distribute queries into blocks. """
    block_size = ceil(len(queries) / num_blocks)
    return [dict(list(queries.items())[i * block_size:(i + 1) * block_size]) for i in range(num_blocks)]

def main():
    parser = argparse.ArgumentParser(description="Align FASTQ reads against a reference genome.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ file containing reads.")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome in FASTA format.")
    parser.add_argument("-k", "--kmer-size", type=int, help="Manually override the k-mer size for indexing and searching.")
    parser.add_argument("-t", "--threshold", type=float, default=0.2, help="Significance threshold for determining significant k-mers based on GC-content.")
    parser.add_argument("-e", "--entropy-threshold", type=float, default=1.5, help="Entropy threshold for determining significant k-mers.")
    parser.add_argument("-n", "--num-threads", type=int, default=4, help="Number of threads to use for processing.")
    
    args = parser.parse_args()
    references = parse_sequences(args.reference, 'fasta')
    queries = parse_sequences(args.input, 'fastq')  # This now returns a dict with (sequence, quality)
    num_threads = args.num_threads  # Number of threads adjusted via command-line
    queries_blocks = distribute_queries(queries, num_threads)
    
    print("@HD\tVN:1.6\tSO:unsorted")
    for ref_id in references:
        print(f"@SQ\tSN:{ref_id}\tLN:{len(references[ref_id])}")

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(align_block, block, references, args.kmer_size, args.threshold, args.entropy_threshold) for block in queries_blocks]
        for future in as_completed(futures):
            for result in future.result():
                print(result)

if __name__ == "__main__":
    main()

