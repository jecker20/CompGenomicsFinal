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
The naive search program reads in synthetic genetic and protein sequences data and stores is it in a dictionary that assigns seqeunce names as keys to the seqeunce within the program. The protein data dictionary is then used to make a dictionary of k-mer indices of each protein sequence with their seqeunce names assigned as keys to their respective k-mer index. The program then takes the genetic data dictionary and goes through each read and translates all 6 possible open reading frames into peptide sequences. These are then searched against the dictionary of protein k-mer indices to find a match using the k-length prefixes as the initial search. If prefixes match, then a length-checking heuristic is applied and if the string lengths do not match, then further comparison is skipped. This is possible under the assumption that all macthes are perfect and complete. If lengths do match, the program continues comparisons with the python "==" operator to check the seqeunces do indeed match. If a mathc is found the seqeunce is removed from the protion k-mer index dictionary. A dicitonary is returned, with open reading frame sequence names (named in-program) as keys to the sequence name in the synthetic protein file. To run the program on terminal, 3 command-line arguments are needed (in order): synthetic DNA data file name, synthetic protein data file name, and an output file name. The program will write an output file with the given name that has all matches on separate lines. k-value is modified direclty in the code as an argument for the main funciton.

# Shared K-mer search

Usage:
```
python less_naive.py [-h] -k <k-mer k-value> -p <proteins file path> -r <DNA read file path> -o <output file path>
```
Example:
```
python .\less_naive.py -k 4 -p .\100-10\100-10-0-protein.txt -r .\100-10\100-10-0-read.txt -o match_output.txt
```

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
```
python k_benchmark.py --kmin 5 --kmax 50 --kstep 2 -l 100 -n 10 -f 10 -r 5 -m shared
```
