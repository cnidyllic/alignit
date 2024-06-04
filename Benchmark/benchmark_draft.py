import os
import subprocess
import time
import matplotlib.pyplot as plt
from memory_profiler import memory_usage

# Define the paths to the input files and scripts
FASTA_FILE = "GCF_000002765.6.fasta"
FASTQ_FILE = "SRR006911.fastq"
ALIGN_IT_SCRIPT = "./align-it/align_it.py"
BWA_MEM_PATH = "bwa"

# Benchmarking function for any alignment tool
def benchmark_tool(command, output_file):
    start_time = time.time()
    peak_memory = max(memory_usage((subprocess.run, (command,), {'stdout': open(output_file, 'w'), 'stderr': subprocess.PIPE})))
    end_time = time.time()

    runtime = end_time - start_time
    return runtime, peak_memory

# Compare align_it with bwa_mem
def compare_alignments():
    # Run align_it
    align_it_command = ["python", ALIGN_IT_SCRIPT, "-i", FASTQ_FILE, "-r", FASTA_FILE]
    align_it_runtime, align_it_memory = benchmark_tool(align_it_command, "align_it_output.sam")

    # Run bwa_mem
    bwa_command = [BWA_MEM_PATH, FASTA_FILE, FASTQ_FILE]
    bwa_runtime, bwa_memory = benchmark_tool(bwa_command, "bwa_mem_output.sam")

    # Plotting the results
    tools = ['align_it', 'bwa_mem']
    runtimes = [align_it_runtime, bwa_runtime]
    memories = [align_it_memory, bwa_memory]

    plt.figure(figsize=(10, 5))
    
    plt.subplot(1, 2, 1)
    plt.bar(tools, runtimes, color=['blue', 'green'])
    plt.xlabel('Tool')
    plt.ylabel('Runtime (seconds)')
    plt.title('Runtime Comparison')

    plt.subplot(1, 2, 2)
    plt.bar(tools, memories, color=['red', 'purple'])
    plt.xlabel('Tool')
    plt.ylabel('Memory Usage (MiB)')
    plt.title('Memory Usage Comparison')

    plt.suptitle('Alignment Tool Benchmark: align_it vs bwa_mem')
    plt.savefig('Benchmark/results_graphs.png')
    plt.show()

def main():
    if not os.path.exists('Benchmark'):
        os.makedirs('Benchmark')
    compare_alignments()

if __name__ == "__main__":
    main()
