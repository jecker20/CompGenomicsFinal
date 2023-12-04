# CompGenomicsFinal

# Data Generator

# Check for Accuracy

# Naive search

# K-mer search

# Benchmarks

### Performance benchmark for read length and number of reads

Usage:
```
python length_benchmark.py -l <read length> -n <reads per file> -f <number of files> -r <repeats per file> -k <k-mer value> -m <naive/shared/random>
```
Example:
```
python length_benchmark.py -l 100 -n 10 -f 10 -r 5 -k 5 -m shared
```

Output is processing time in seconds, one line per file, averaged across all repeats of that file.
