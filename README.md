# align_it (CSE185 Project Demo)
This is a demonstration project for CSE 185. It is a Python-based tool designed for aligning reads to the Influenza A reference genome. This tool is optimized for handling complex genomic sequences with high AT/low GC content and repetitive elements specifically.
# Features
- Dynamically adjusts k-mer size based on GC content and entropy
- Builds a hash-based index for k-mers from the reference sequence
- Aligns FASTQ reads to the reference genome and outputs results in SAM format
- Multithreading support for faster processing
# Installation Instructions
Clone the repository:
```
git clone https://github.com/cnidyllic/alignit
```
Change into directory and install tool so it can be used from command line:
```
cd alignit
```
```
python setup.py install
```
Note: if you do not have root access, you can run the commands above with additional options to install locally: 
```
python setup.py install –-user
```
If the install was successful, typing ```python align-it/align_it.py --help``` should show a useful message. 

# Basic usage: 
The basic usage of ```align_it``` is:
```
python align-it/align_it.py -i <input_FASTQ_file> -r <reference_genome>
```
To run ```align_it``` on a small test example:
```
python align-it/align_it.py -i example-files/test_queries.fastq -r example-files/test_reference.fa 
```
Complete usage instructions: 

# align_it options
The required inputs to ```align_it``` is a FASTQ file containing the reads and the reference genome fasta file for alignment. Users may additionally specify the options below:

-i, --input, required, "Input FASTQ file containing reads."

-r, --reference, required, "Reference genome in FASTA format."

-k, --kmer-size, type=int, "Manually override the k-mer size for indexing and searching."

-t, --threshold, type=float, default=0.2, "Significance threshold for determining significant k-mers based on GC-content."

-e, --entropy-threshold, type=float, default=1.5, "Entropy threshold for determining significant k-mers."

-n, --num-threads, type=int, default=4, "Number of threads to use for processing, or number of logical cores."

# Benchmarking:
Install necessary tools first:
```pip install memory_profiler```

# File format
The output file format is the same as the bwa mem method, a SAM file, printed to the standard output: See: python align-it/align_it.py -i example-files/test_queries.fastq -r example-files/test_reference.fa 

# Citations
This project references the following works:

Altschul, S.F., Gish, W., Miller, W., Myers, E.W., & Lipman, D.J. (1990). Basic local alignment search tool. Journal of Molecular Biology, 215(3), 403-410.
- Hash-based indexing inspired from this study.

Langmead, B., Trapnell, C., Pop, M., & Salzberg, S.L. (2009). Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biology, 10(3), R25.
- Mentions the use of k-mer sizes to align sequences.

Compeau, P.E.C., Pevzner, P.A., & Tesler, G. (2011). How to apply de Bruijn graphs to genome assembly. Nature Biotechnology, 29(11), 987-991.
- Detailing k-mer size use specifically.

# License
This project is licensed under the MIT License. See the LICENSE file for details.

# Contributors
This repository was generated by: Corey Nguyen, with inspiration from https://github.com/gymreklab/cse185-demo-project#readme, and work of Jenny Nguyen and Tiffany Zhang. 
Please submit a pull request with any corrections and suggestions.

