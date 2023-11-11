#AOA
# Modified code form solutions provided for hw2
import sys
from collections import defaultdict

def parse_syntheticDNA(fn):
    reads = []
    names = []
    with open(fn, 'rt') as fh:
        while True:
            name = fh.readline().strip()
            if len(name) == 0:
                break
            names.append(name)
            st = fh.readline().strip()
            reads.append(st)
    return reads, names

def parse_syntheticProt(fn):
    reads = []
    names = []
    with open(fn, 'rt') as fh:
        while True:
            name = fh.readline().strip()
            if len(name) == 0:
                break
            if name[0] != '>': # check that its not junk
                break
            name = name[1:]
            names.append(name)
            st = fh.readline().strip()
            reads.append(st)
    return reads, names

def make_index(text, k):
    index = defaultdict(list)
    for i in range(len(text) - k + 1):
        substr = text[i:i+k]
        index[substr].append(i)
    return index


def make_dataset_index(reads, names, k):
    if len(reads) != len(names):
        raise Exception("#names != #reads")
    index = {}
    for i in range(len(reads)): # this has to be for data in set
        index[names[i]] = make_index(reads[i], k)
    return index

# CURRENTL NON-FUNCTIONAL
#EXACT MATCHING: filters out size differences
def exact_matches(peptides, indexDict, pepts, k):
    matches = defaultdict(list)
    for p in peptides:
        for n in indexDict.keys():
            if p[:k] in indexDict[n]:
                for o in indexDict[n][p[:k]]:
                    if genome[o:o+len(p)] == p:
                        matches[p].append(o)
                        matches[p] = list(set(matches[p]))
    return matches



def revComp(seq):
    str = ""
    tmp = seq[::-1]

    for n in tmp:
        match n:
            case 'A':
                str += 'T'
            case 'T':
                str += 'A'
            case 'C':
                str += 'G'
            case 'G':
                str += 'C'
    return str

def make_orfs(seq):
    revSeq = revComp(seq)
    orfs = [seq, seq[1:], seq[2:], revSeq, revSeq[1:], revSeq[2:]]
    return orfs

def transcript(seq):
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
    seq = transcript(seq)
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)] # separate codons
    pep = ""
    for c in codons:
        if len(c) < 3:
            break
        #print(c)
        pep += codon_match(c)
    #print(pep)
    
    return pep

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



if __name__ == "__main__":
    dnaFile = sys.argv[1]
    protFile = sys.argv[2]
    infile = sys.argv[3]
    outFile = sys.argv[4]

    readK = 4
    protK = 4

    reads, dNames = parse_syntheticDNA(dnaFile)

    pepts, pNames = parse_syntheticProt(protFile)

    readIndex = make_dataset_index(reads, dNames, readK)
    protIndex = make_dataset_index(pepts, pNames, protK)

    with open(infile, 'r') as fh:
        p = fh.readline().strip()

    orfs = make_orfs(p)
    peptides = []

    for frame in orfs:
        peptides.append(translate(frame))

    #print(translate(reads[4]))

    print(orfs)
    print(peptides)

    