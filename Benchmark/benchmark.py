import time
import matplotlib.pyplot as plt
import numpy as np
from memory_profiler import memory_usage
from src.main import main as align_main
from metrics import calculate_accuracy

def run_benchmark(input_path, reference_path):
    # Measure runtime and memory usage
    start_time = time.time()
    memory_usage_profile = memory_usage((align_main, (input_path, reference_path)), interval=0.1, include_children=True)
    end_time = time.time()

    runtime = end_time - start_time
    peak_memory = max(memory_usage_profile)

    # Generate output path based on input for automated accuracy measurement
    output_path = input_path.replace('.fastq', '_aligned.sam')

    # Calculate accuracy
    accuracy = calculate_accuracy(output_path, 'expected_output_path')  # Placeholder

    return runtime, peak_memory, accuracy, memory_usage_profile

def plot_results(times, memories, accuracies):
    # Time and memory plot
    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('Run #')
    ax1.set_ylabel('Runtime (s)', color=color)
    ax1.plot(times, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:blue'
    ax2.set_ylabel('Memory Usage (MiB)', color=color)  # we already handled the x-label with ax1
    ax2.plot(memories, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.title('Runtime and Memory Usage over Runs')
    plt.show()

    # Accuracy plot
    plt.figure()
    plt.plot(accuracies, label='Accuracy', marker='o')
    plt.xlabel('Run #')
    plt.ylabel('Accuracy (%)')
    plt.title('Alignment Accuracy over Runs')
    plt.legend()
    plt.show()

def main():
    input_paths = ["data/query_sequences/sample1.fastq", "data/query_sequences/sample2.fastq"]
    reference_path = "data/reference_genomes/sample.fasta"
    times, memories, accuracies = [], [], []

    for input_path in input_paths:
        runtime, peak_memory, accuracy, _ = run_benchmark(input_path, reference_path)
        times.append(runtime)
        memories.append(peak_memory)
        accuracies.append(accuracy * 100)  # convert to percentage

    plot_results(times, memories, accuracies)

if __name__ == "__main__":
    main()
