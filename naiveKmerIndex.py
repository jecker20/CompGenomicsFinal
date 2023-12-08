# NAIVE KMER INDEX SEARCHING PROGRAM 

import sys
from collections import defaultdict


def parse_syntheticDNA(fn):
    """ Method to parse through a synthetic DNA reads data file and
    store into a dictionary that assigns read names as keys to 
    DNA sequence reads

    Args:
        fn (string): a string filename of the synthetic DNA reads data file

    Returns:
        dct (dictionary): A dictionary with the read names as keys
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


def parse_syntheticProt(fn):
    """ Method to parse through a synthetic protein sequences data file and
    store into a dictionary that assigns protein names as keys to 
    peptide sequences

    Args:
        fn (string): filename of the synthetic protein reads data file

    Returns:
        dct (dictionary): A dictionary with the peptide names as keys
            to their respective sequences
    """
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


def reverse(seq):
    """ Helper function that reverses a string, used to find
    all open reading frames

    Args:
        seq (string): a DNA seqeunce of interest

    Returns:
        str (string): the reverse of the inputted DNA sequence
    """
    str = seq[::-1]

    return str

def make_orfs(seq):
    """ Helper function that returns list of all 6 possible ORFs
    of a DNA sequence

    Args:
        seq (str): A DNA sequence of interest

    Returns:
        orfs (list): a list of strings of all 6 possible ORFs
    """
    revSeq = reverse(seq)
    orfs = [seq, seq[1: len(seq) - 2], seq[2: len(seq) - 1], revSeq, revSeq[1: len(seq) - 2], revSeq[2: len(seq) - 1]]
    return orfs

def transcribe(seq):
    """ Helper function that returns a the complement RNA sequence to a DNA
    sequence of interest

    Args:
        seq (string): DNA sequence of interest

    Returns:
        rna (string): compliment RNA sequence
    """
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

def translate(seq):
    """ Helper function that returns a peptide sequence of a given
    RNA sequence of interest

    Args:
        seq (string): RNA sequence of interest

    Returns:
        pep (string): peptide sequence derived from inputted RNA sequence
    """
    seq = transcribe(seq)
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)] # separate codons
    pep = ""
    for c in codons:
        if len(c) < 3:
            break
        pep += codon_match(c)
    
    return pep

def codon_match(codon):
    """ helper fucntions that takes a codon and returns
    its appropiate amino acid translation

    Args:
        codon (string): a 3-length string of an RNA codon

    Returns:
        amino (char): a character that represents the 1-letter abbrevation
        of the translated amino
    """
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
        print(len(protIndexDict.keys()))
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

    
