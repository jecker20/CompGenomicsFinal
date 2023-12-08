# CompGenomicsFinal

# Data Generator
Used to generate data files for analysis. 
Creates two files
1) Read file 
   A) Contains genomic reads
   B) Reads labelled as >read#
2) Protein file
   A) Contains translated protein assocaited with generated reads
   B) Contains junk protein
   C) Proteins connected to reads labelled >read#-# where first # is assocaited read and second # is assocaited reading frame
   D) Junk proteins labelled >junkread#
Files naming scheme 'length of read'-'number of reads'-'file number'-'read/protein'.txt
Variables for adjusting data generation
1) LENGTH = length of read
2) numberOfFiles = number of files to generate
3) NUMOFREADS = number of reads to generate
Amount of junk and kept protein is set to a random variable but this can be changed to set integers if desired (junk is always greater than 1 which provides files with more junk than matches)


# Check for Accuracy
This program confirms that a matching program provided a correct output
Output: The first line of output indicates if the read-match pairs are all correct. The second line of output indciates if the number of matches equals the number of expected matches based on the provided protein file
Input: Variables matchedFile and proteinFile are modifiable vairbales that must match (protein file must be the file that match was run against). The matchedFile must be in a new line space list of "read,match"

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

### Performance benchmark to compare k-mer k values

Usage
```
python k_benchmark.py --kmin <value> --kmax <value> --kstep <value> -l <read length> -n <reads per file> -f <number of files> -r <repeats per file> -k <k-mer value> -m <naive/shared/random>
```

Example:

