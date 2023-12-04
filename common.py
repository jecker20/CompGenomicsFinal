from collections import defaultdict

# Parses through SYNTHETIC DNA reads file and stores reads
# into a dict w/ names as keys
def parse_syntheticDNA(fn):
    dct = {}
    with open(fn, 'rt') as fh:
        while True:
            name = fh.readline().strip()
            if len(name) == 0:
                break
            st = fh.readline().strip()
            dct[name] = st
    return dct 

# Parses through SYNTHETIC TRANSLATED PETIDES file and stores reads
# into a dict w/ names as keys
def parse_syntheticProt(fn):
    dct = {}
    with open(fn, 'rt') as fh:
        while True:
            name = fh.readline().strip()
            if len(name) == 0:
                break
            if name[0] != '>': # check that its not junk (synthetic property)
                break
            name = name[1:]
            st = fh.readline().strip()
            dct[name] = st
    return dct

# makes a kmer index - derived from solution in homework 2
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

# Helper function: returns list of 6 possible ORFs of a DNA sequence
def make_orfs(seq):
    revSeq = seq[::-1]
    orfs = [seq, seq[1:-2], seq[2:-1], revSeq, revSeq[1:-2], revSeq[2:-1]]
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
