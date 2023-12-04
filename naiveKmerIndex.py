# NAIVE KMER INDEX SEARCHING PROGRAM 

import sys
from collections import defaultdict

from common import parse_syntheticProt, parse_syntheticDNA, make_orfs, translate

def make_index(text, k):
    index = defaultdict(list)
    for i in range(len(text) - k + 1):
        substr = text[i:i+k]
        index[substr].append(i)
    return index

# makes a dict of kmer indexs for multiple reads/proteins
def make_dataset_index(reads, k):
    indexDict = {}
    for n in reads.keys():
        indexDict[n] = make_index(reads[n], k)
    return indexDict

# MAIN NAIVE LOOKUP FUNCTION: takes all 6 possible ORF peptide sequences
# and finds if their is matches in synProtein data kmerIndex dictionary
# DISCLAIMER: does NOT interact with DNA read kmer indexs of any sort.
def exact_matches(read, peptides, protIndexDict, synProts, k):
    matches = defaultdict()
    for i in range(len(peptides)):
        for n in protIndexDict.keys():
            currIdx = protIndexDict[n]
            if peptides[i][:k] in currIdx:
                name = n
                if len(peptides[i]) != len(synProts[name]):
                    pass
                else:
                    if peptides[i] == synProts[name]:
                        matches[str(read) + '-' + str(i + 1)] = n
    return matches

# Main function called, seperated for benchmarking
def main(kVal, protFile, dnaFile, output_file=None, repeats=1):
    # Read Input
    synReads = parse_syntheticDNA(dnaFile)
    synProts = parse_syntheticProt(protFile)

    # Make Dateset Index
    protIndex = make_dataset_index(synProts, kVal)
    #readIndex = make_dataset_index(synReads, readK) # more relevant to "shared-kmer" searching

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
    main(5, protFile, dnaFile, outFile)

