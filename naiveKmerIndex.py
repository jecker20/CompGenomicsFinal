# NAIVE KMER INDEX SEARCHING PROGRAM 

import sys
from collections import defaultdict


def parse_syntheticDNA(fn):
    """_summary_

    Args:
        fn (function): _description_

    Returns:
        dct : A dictionary with the read names as keys
            to their respective sequences
    """
    dct = {}
    with open(fn, 'rt') as fh:
        while True:
            name = fh.readline().strip()
            if len(name) == 0:
                break
            st = fh.readline().strip()
            dct[name] = st
    return dct 

# Parses through SYNTHETIC PEPTIDE SEQUENCES file and stores reads
# into a dict w/ names as keys
def parse_syntheticProt(fn):
    dct = {}
    with open(fn, 'rt') as fh:
        while True:
            name = fh.readline().strip()
            if len(name) == 0:
                break
            if name[0] != '>': # check its data (synthetic property)
                break
            name = name[1:] # remove arrowhead
            st = fh.readline().strip()
            dct[name] = st
    return dct

""" 
"""
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

# Helper function: returns reverse compliment of a DNA sequence
def reverse(seq):
    str = seq[::-1]

    return str

# Helper function: returns list of 6 possible ORFs of a DNA sequence
def make_orfs(seq):
    revSeq = reverse(seq)
    orfs = [seq, seq[1: len(seq) - 2], seq[2: len(seq) - 1], revSeq, revSeq[1: len(seq) - 2], revSeq[2: len(seq) - 1]]
    return orfs

# Helper function: returns transcript of a DNA sequence into RNA
# to aid translation
def transcribe(seq):
    rna = ""
    for b in seq:
        match b:
            case 'A':
                rna += 'U'
            case 'C':
                rna += 'G'
            case 'G':
                rna += 'C'
            case 'T':
                rna += 'A'
    
    return rna

# Helper function: returns peptide sequence of a given RNA sequence
def translate(seq):
    seq = transcribe(seq)
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)] # separate codons
    pep = ""
    for c in codons:
        if len(c) < 3:
            break
        pep += codon_match(c)
    
    return pep

# Helper function: translates matches a codon to appropiate amino
def codon_match(codon):
    amino = ''
    match codon[0]: # read first base
        case 'A':
            match codon[1]: # read second base
                case 'A':
                    match codon[2]: # read third base
                        case 'A' | 'G':
                            amino = 'K' # Lys (K)
                        case 'C' | 'U':
                            amino = 'N' # Asn (N)
                case 'C': # cases where not all bases need to be read
                    amino = 'T'   # Thr (T)
                case 'G':
                    match codon[2]:
                        case 'A' | 'G':
                            amino = 'R' # Arg (R)
                        case 'C' | 'U':
                            amino = 'S' # Ser (S)
                case 'U':
                    match codon[2]:
                        case 'G':
                            amino = 'M' # Met (M)
                        case 'A' | 'C' | 'U':
                            amino = 'I' # Ile (I)
        
        case 'C':
            match codon[1]:
                case 'A':
                    match codon[2]:
                        case 'A' | 'G':
                            amino = 'Q' # Gln (Q)
                        case 'C' | 'U':
                            amino = 'H' # His (H)
                case 'C':
                    amino = 'P' # Pro (P)
                case 'G':
                    amino = 'R' # Arg (R)
                case 'U':
                    amino = 'L' # Leu (L)
        
        case 'G':
            match codon[1]:
                case 'A':
                    match codon[2]:
                        case 'A' | 'G':
                            amino = 'E' # Glu (E)
                        case 'C' | 'U':
                            amino = 'D' # Asp (D)
                case 'C':
                    amino = 'A' # Ala (A)
                case 'G':
                    amino = 'G' # Gly (G)
                case 'U':
                    amino = 'V' # Val (V)
        
        case 'U':
            match codon[1]:
                case 'A':
                    match codon[2]:
                        case 'A' | 'G':
                            amino = '-' # STOP
                        case 'C' | 'U':
                            amino = 'Y' # Tyr (Y)
                case 'C':
                    amino = 'S' # Ser (S)
                case 'G':
                    match codon[2]:
                        case 'A':
                            amino = '-' # STOP
                        case 'C' | 'U':
                            amino = 'C' # Cys (C)
                        case 'G':
                            amino = 'W' # Trp (W)
                case 'U':
                    match codon[2]:
                        case 'A' | 'G':
                            amino = 'L' # Leu (L)
                        case 'C' | 'U':
                            amino = 'F' # Phe (F)
        
        case _:
            pass
    
    return amino

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
def main(kVal, protFile, dnaFile, outFile):
    # Read Input
    synReads = parse_syntheticDNA(dnaFile)
    synProts = parse_syntheticProt(protFile)

    # Make Dateset Index
    protIndex = make_dataset_index(synProts, kVal)
    #readIndex = make_dataset_index(synReads, readK) # more relevant to "shared-kmer" searching
    
    # Output
    with open(outFile, 'w') as outh:
        for k in synReads.keys():
            orfs = make_orfs(synReads[k])
            peptides = []
            for frame in orfs:
                peptides.append(translate(frame))
            mats = exact_matches(k, peptides, protIndex, synProts, kVal)
            for k in mats.keys():
                outh.write(str(k) + ',' + mats[k] + '\n')


if __name__ == "__main__":
    # Input
    dnaFile = sys.argv[1]
    protFile = sys.argv[2]
    outFile = sys.argv[3]

    # Call Main Function
    main(5, protFile, dnaFile, outFile)

