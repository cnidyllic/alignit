import os
import sys
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from math import log2, ceil
import time
from memory_profiler import memory_usage


# Performs alignment of sequences based on kmer
class align_it:
    # Initialize an instance of the align_it 
    def __init__(self, reference_seq, k_override=None, significance_threshold=0.2, entropy_threshold=1.5):
        self.reference_seq = reference_seq # store the reference sequence
        self.significance_threshold = significance_threshold
        self.entropy_threshold = entropy_threshold 
        self.k = k_override or 20 # kmer size, default to 20 unless overridden in option
        self.index = self.build_index(self.reference_seq, self.k) # builds an index of kmers from the reference

    # Calculate Shannon entropy for a kmer
    def calculate_entropy(self, kmer):
        freq = {x: kmer.count(x) / len(kmer) for x in set(kmer)} # calculate freq of each base in the kmer
        entropy = -sum(p * log2(p) for p in freq.values() if p > 0) # computes shannon entropy
        return entropy

    # adjusts kmer size based on GC content
    def adjust_kmer_size(self, sequence):
        if self.k is not None:
            return self.k 
        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence) # calculate GC content
        if gc_content < 0.4:
            return max(20, len(sequence) // 100)
        elif gc_content > 0.6:
            return max(15, len(sequence) // 150)
        return 20

    # Hash-based index for kmers from the reference sequence
    def build_index(self, sequence, k):
        index = {}
        step = max(1, k // 4) # determine the step size for sliding the kmer to reduce overlaps in AT-rich seq
        # Iterate over reference sequnce to create kmers and stores positions in index
        for i in range(0, len(sequence) - k + 1, step):
            kmer = sequence[i:i+k]
            if self.is_significant_kmer(kmer): 
                index.setdefault(kmer, []).append(i)
        return index

    # Signifance of kmer based on GC content and entropy
    def is_significant_kmer(self, kmer):
        gc_content_significant = (kmer.count('G') + kmer.count('C')) / len(kmer) > self.significance_threshold
        entropy_significant = self.calculate_entropy(kmer) > self.entropy_threshold
        return gc_content_significant and entropy_significant

    # Searches for the read in the reference sequence using the kmer index
    def search(self, query):
        k = self.adjust_kmer_size(query) # adjust kmer size based on raw reads
        self.index = self.build_index(self.reference_seq, k) # rebuild the index with adjusted kmer size
        matches = set()
        # Search for each kmer of raw reads in the index
        for i in range(len(query) - k + 1):
            kmer = query[i:i+k]
            if kmer in self.index:
                for pos in self.index[kmer]:
                    end_pos = pos + len(query)
                    if end_pos <= len(self.reference_seq) and query == self.reference_seq[pos:end_pos]:
                        matches.add(pos)
        return sorted(list(matches))

# Parsing from a correctly formatted file, FASTA or FASTQ only
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

def quality_string_to_scores(quality_string):
    return [ord(char) - 33 for char in quality_string]

def summarize_quality(quality_string):
    scores = [ord(char) - 33 for char in quality_string]
    average_score = sum(scores) / len(scores)
    return f"Average Quality: {average_score:.2f}"
    
# align a block of reads to a set of reference sequences
def align_block(queries_block, references, k_override, significance_threshold, entropy_threshold):
    results = []
    aligned_count = 0  # To count the number of successfully aligned reads
    # iterate over each query and reference to align
    for reference_id, reference_seq in references.items():
        aligner = align_it(reference_seq, k_override, significance_threshold, entropy_threshold)
        for query_id, query_info in queries_block.items():
            query_seq, query_qual = query_info
            match_positions = aligner.search(query_seq)
            quality_summary = summarize_quality(query_qual)
            if match_positions:
                aligned_count += 1
                for pos in match_positions:
                    results.append(f"{query_id}\t0\t{reference_id}\t{pos+1}\t255\t{len(query_seq)}M\t*\t0\t0\t{query_seq}\t{quality_summary}")
            else:
                results.append(f"{query_id}\t4\t*\t0\t0\t*\t*\t0\t0\t{query_seq}\t{quality_summary}")
    return results, aligned_count

# distribute queries into block for processing (threading)
def distribute_queries(queries, num_blocks):
    block_size = ceil(len(queries) / num_blocks) # calculate size of each block (evenly)
    # split the reads into the specified number of blocks
    return [dict(list(queries.items())[i * block_size:(i + 1) * block_size]) for i in range(num_blocks)]

def main():
    parser = argparse.ArgumentParser(description="Align FASTQ reads against a reference genome.")
    # command-line arguments for the script
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ file containing reads.")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome in FASTA format.")
    parser.add_argument("-k", "--kmer-size", type=int, help="Manually override the k-mer size for indexing and searching.")
    parser.add_argument("-t", "--threshold", type=float, default=0.2, help="Significance threshold for determining significant k-mers based on GC-content.")
    parser.add_argument("-e", "--entropy-threshold", type=float, default=1.5, help="Entropy threshold for determining significant k-mers.")
    parser.add_argument("-n", "--num-threads", type=int, default=4, help="Number of threads to use for processing, or number of logical cores.")
    
    args = parser.parse_args() # parse command line arguments    
    
    # Start measuring time and memory
    start_time = time.time()
    mem_usage_before = memory_usage(max_usage=True)

    references = parse_sequences(args.reference, 'fasta') # parse reference genome
    queries = parse_sequences(args.input, 'fastq') # parse raw reads
    num_threads = args.num_threads # number of threads for parallel processing
    queries_blocks = distribute_queries(queries, num_threads) # distribute reads into blocks
    total_aligned = 0
    total_reads = sum(len(block) for block in queries_blocks)

    output_file = open("output.sam", "w")
    output_file.write("@HD\tVN:1.6\tSO:unsorted\n")
    for ref_id in references:
        output_file.write(f"@SQ\tSN:{ref_id}\tLN:{len(references[ref_id])}\n")

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(align_block, block, references, args.kmer_size, args.threshold, args.entropy_threshold) for block in queries_blocks]
        results, aligned_count = future.result()
        total_aligned += aligned_count
        for result in results:
            output_file.write(result + "\n")

    output_file.close()

    # Calculate and print execution time and memory usage
    end_time = time.time()
    runtime = end_time - start_time
    mem_usage_after = memory_usage(max_usage=True)
    peak_memory_usage = mem_usage_after - mem_usage_before

    print(f"Runtime: {runtime:.2f} seconds")
    print(f"Peak Memory Usage: {peak_memory_usage:.2f} MiB")
    print(f"Total Reads: {total_reads}, Reads Aligned: {total_aligned}")
    print(f"Alignment Success Rate: {alignment_rate:.2%}")

if __name__ == "__main__":
    main() # if python script is executed directly
