# NAIVE KMER INDEX SEARCHING PROGRAM 

import sys
from collections import defaultdict

from common import parse_syntheticDNA, parse_syntheticProt, make_orfs, translate

def make_index(text, k):
    """ Method to create k-mer index of a given string, inspired by
    method provided in homework 2 answers

    Args:
        text (string): string that is being turned into a k-mer index
        k (int): value for length of k-mers

    Raises:
        Exception: do not make index if k is greater than string length

    Returns:
        index (dictionary): k-mer index stored into a dicitonary
    """
    if k >= len(text):
        raise Exception("k value is larger than sequence length")
    index = defaultdict(list)
    for i in range(len(text) - k + 1):
        substr = text[i:i+k]
        index[substr].append(i)
    return index

def make_dataset_index(reads, k):
    """ Method to make a dictionary that assigns a sequence name as a key to
    their respecitve k-mer index, used to make k-mer indices for all data
    in a set

    Args:
        reads (dictionary): dictionary with names as keys to sequences
        k (int): value for k-mer length

    Returns:
        indexDict (dictionary): a dictionary that assigns a read name as a key 
        to its respecitve seqeunce k-mer index
    """
    indexDict = {}
    for n in reads.keys():
        indexDict[n] = make_index(reads[n], k)
    return indexDict

def exact_matches(read, peptides, protIndexDict, synProts, k):
    """ MAIN NAIVE LOOKUP FUNCTION: takes all 6 possible ORF peptide sequences
    and finds if their is matches in synthetic protein data kmerIndex dictionary
    DISCLAIMER: does NOT interact with DNA read kmer indexs of any sort.
    NOTE: inspired by homework 2

    Args:
        read (string): DNA sequence NAME of interest for attempting to find matches in database
        peptides (list): list of all possible peptide sequences from 6 possibel ORFs
        protIndexDict (dicitonary): dictionary that assigns a peptide seqeunce's name as a key
        to its respective k-mer index
        synProts (dicitionary): dictionary that assigns a peptide seqeunce's name as a key
        to its respective whole peptide sequence
        k (int): value for length of k-mers

    Returns:
        matches (dictionary): dictioanry that assigns peptide names as keys to their
        respective matches
    """
    matches = defaultdict()
    for i in range(len(peptides)):
        for n in protIndexDict.copy().keys(): # copy allows to remove matches from dict once found
            currIdx = protIndexDict[n]
            if peptides[i][:k] in currIdx:
                name = n
                if len(peptides[i]) != len(synProts[name]):
                    pass
                else:
                    if peptides[i] == synProts[name]:
                        matches[str(read) + '-' + str(i + 1)] = n # naming convention
                        protIndexDict.pop(n) # remove from future searching if match found
    return matches

def main(kVal, protFile, dnaFile, output_file=None, repeats=1):
    """ Main funtion that runs the naive searching program

    Args:
        kVal (int): value for k-mer lengths
        protFile (string): file name of synthetic peptide sequences
        dnaFile (string): file name of synthetic DNA reads
        output_file (string, optional): file name for written output file. Defaults to None.
        repeats (int, optional): integer value for number of times to run search method; used for
        benchmarking trials. Defaults to 1.

    Output:
        Writes output to the output file name provided
    """

    # Read Input
    synReads = parse_syntheticDNA(dnaFile)
    synProts = parse_syntheticProt(protFile)

    # Make Dateset Index
    protIndex = make_dataset_index(synProts, kVal)
    # no DNA index as it is more relevant to "shared-kmer" searching

    for _ in range(repeats):
        matches = {}

        for k in synReads.keys():
            orfs = make_orfs(synReads[k])
            peptides = []
            for frame in orfs:
                peptides.append(translate(frame))
            mats = exact_matches(k, peptides, protIndex, synProts, kVal)
            for k in mats.keys():
                matches[k] = mats[k]

    # Output
    if output_file is not None:
        with open(output_file, 'w') as outh:
           for k in matches.keys():
                    outh.write(str(k) + ',' + matches[k] + '\n')


if __name__ == "__main__":
    # Input
    dnaFile = sys.argv[1]
    protFile = sys.argv[2]
    outFile = sys.argv[3]
    # Call Main Function
    main(100, protFile, dnaFile, outFile)
