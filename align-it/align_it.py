import os
import sys
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from math import log2, ceil
import time
from memory_profiler import memory_usage


# Performs alignment of sequences based on kmer
class align_it:
    # Constructor: Initializes the alignment tool with reference sequence and parameters.
    def __init__(self, reference_seq, k_override=None, significance_threshold=0.4, entropy_threshold=1.5):
        self.reference_seq = reference_seq 
        self.significance_threshold = significance_threshold 
        self.entropy_threshold = entropy_threshold 
        self.k = k_override or 20 # kmer size, default to 20 unless overridden in option
        self.index = self.build_index(self.reference_seq, self.k) # builds an index of kmers from the reference

     # Calculates the Shannon entropy of a k-mer to evaluate its complexity.
    def calculate_entropy(self, kmer):
        freq = {x: kmer.count(x) / len(kmer) for x in set(kmer)} # calculate freq of each base in the kmer
        entropy = -sum(p * log2(p) for p in freq.values() if p > 0) # computes shannon entropy
        return entropy

    # Adjusts kmer size based on GC content
    def adjust_kmer_size(self, sequence):
        if self.k is not None:
            return self.k # Returns overridden k-mer size if set.
        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence) # calculate GC content
        # Adjusts k-mer size based on GC content thresholds.
        if gc_content < 0.4:
            return max(15, len(sequence) // 200)
        elif gc_content > 0.6:
            return max(25, len(sequence) // 150)
        return 20

    # Hash-based index for kmers from the reference sequence
    def build_index(self, sequence, k):
        index = {}
        step = max(1, k // 4) # Step size for reducing overlaps in AT-rich sequences.
        # Iterate over reference sequnce to create kmers and stores positions in index
        for i in range(0, len(sequence) - k + 1, step):
            kmer = sequence[i:i+k]
            if self.is_significant_kmer(kmer): 
                index.setdefault(kmer, []).append(i) # Adds significant k-mers to the index.
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
        sequences = {} # Initialize an empty dictionary to store sequence data.
        with open(file_path, 'r') as file:
            identifier, sequence, quality = '', '', ''
            for line in file: # Iterate over each line in the file.
                line = line.strip() # Strip whitespace from the beginning and end of the line.
                if file_type == 'fasta':
                    if line.startswith('>'):
                        if identifier: # If there's a current identifier, store the previous sequence.
                            sequences[identifier] = sequence 
                        identifier = line[1:] # Set the new identifier, removing the '>'.
                        sequence = '' # Reset sequence for the new entry.
                    else: # If it's not a header, it's part of the sequence.
                        sequence += line
                elif file_type == 'fastq':
                    if line.startswith('@') and not identifier:
                        identifier = line[1:] # Set the new identifier, removing the '@'.
                    elif line.startswith('+') and identifier:
                        continue
                    elif not sequence and not quality:
                        sequence += line # This line must be sequence data, add it.
                    elif sequence and not quality:
                        quality += line # This line must be quality scores, add it.
                    if quality and len(quality) >= len(sequence):
                        sequences[identifier] = (sequence, quality)
                        identifier, sequence, quality = '', '', '' # Reset for the next record.
            if identifier:
                if file_type == 'fasta':
                    sequences[identifier] = sequence # Store the last read sequence.
                else: # If it's FASTQ.
                    sequences[identifier] = (sequence, quality)  # Store the last read sequence and quality.
        return sequences # Return the dictionary containing all sequences.
    except Exception as e:
        sys.stderr.write(f"ERROR: Failed to parse {file_type.upper()} file {file_path}: {str(e)}\n")
        sys.exit(1)

# Converts quality string to numerical scores.
def quality_string_to_scores(quality_string):
    return [ord(char) - 33 for char in quality_string] 

# Calculates and returns the average quality score.
def summarize_quality(quality_string):
    scores = [ord(char) - 33 for char in quality_string]
    average_score = sum(scores) / len(scores)
    return f"Average Quality: {average_score:.2f}"
    
# Align a block of reads to a set of reference sequences using multi-threading.
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

# Distribute queries into block for processing (threading)
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
    
    # Start measuring time and memory for parsing reference sequences.
    start_time_parsing = time.time()
    mem_usage_before = memory_usage(max_usage=True)

    references = parse_sequences(args.reference, 'fasta') # Parse reference genome
    parsing_duration = time.time() - start_time_parsing
    print(f"Parsing Reference Duration: {parsing_duration:.2f} seconds")
    
    # Start measuring time for parsing query sequences.
    start_time_query_parsing = time.time()
    queries = parse_sequences(args.input, 'fastq') # parse raw reads
    query_parsing_duration = time.time() - start_time_query_parsing
    print(f"Parsing Queries Duration: {query_parsing_duration:.2f} seconds")
    
    num_threads = args.num_threads # number of threads for parallel processing
    queries_blocks = distribute_queries(queries, num_threads) # distribute reads into blocks
    
    start_time_alignment = time.time()
    total_aligned = 0
    total_reads = sum(len(block) for block in queries_blocks)

    # Open output file for writing alignments.
    output_file = open("output.sam", "w") 
    output_file.write("@HD\tVN:1.6\tSO:unsorted\n")
    for ref_id in references:
        output_file.write(f"@SQ\tSN:{ref_id}\tLN:{len(references[ref_id])}\n")

    # Create an instance of ThreadPoolExecutor as 'executor'
    # 'max_workers=num_threads' specifies the number of threads to use in the pool
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
         # List comprehension to create future tasks:
        # 'executor.submit' schedules the 'align_block' function to be executed by the pool
        # Each 'block' of queries and the reference genome data are passed as arguments
        # This allows each thread in the pool to process part of the data independently
        futures = [executor.submit(align_block, block, references, args.kmer_size, args.threshold, args.entropy_threshold) for block in queries_blocks]
        
        # Iterate over the futures as they complete (as_completed(futures)):
        # This ensures that processing proceeds as threads finish, regardless of order
        for future in as_completed(futures):
            # Results include the aligned sequences and count of alignments per block
            results, aligned_count = future.result() 
            # Total aligned reads are accumulated from each block processed
            total_aligned += aligned_count
            # Append each result from the block to the output file
            # Each result is a formatted string representing a line in a SAM file
            for result in results:
                output_file.write(result + "\n")

    output_file.close() # Close the output file after writing all results.
    
    alignment_duration = time.time() - start_time_alignment

    print(f"Alignment Processing Duration: {alignment_duration:.2f} seconds")

    # Calculate and print execution time and memory usage
    runtime = time.time() - start_time_parsing
    print(f"Total Runtime: {runtime:.2f} seconds")
    
    mem_usage_after = memory_usage(max_usage=True)
    peak_memory_usage = mem_usage_after - mem_usage_before
    print(f"Peak Memory Usage: {peak_memory_usage:.2f} MiB")

if __name__ == "__main__":
    main() # Start the main function if this script is run as the main program.
