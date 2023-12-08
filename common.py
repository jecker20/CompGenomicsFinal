"""Common helper/utility functions used across code base"""

def parse_syntheticDNA(fn):
    """Parse synthetic DNA reads from file into dict indexed by read name.

    Args:
        fn (str): read file name

    Returns:
        Dict[str, str]: Dict with read name as key, DNA sequence as value
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
    """Parse synthetic proteins from file into dict indexed by protein name.

    Args:
        fn (str): protein file name

    Returns:
        Dict[str, str]: Dict with protein name as key, amino acid sequence as value
    """
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

def make_orfs(seq):
    """Given nucleo sequence, produces 6 possible open reading frames

    Args:
        seq (str): input DNA/RNA sequence, length must be multiple of three!

    Returns:
        List[str]: list of open reading frames
    """
    revSeq = seq[::-1]
    orfs = [seq, seq[1:-2], seq[2:-1], revSeq, revSeq[1:-2], revSeq[2:-1]]
    return orfs

def transcribe(seq):
    """Transcribes DNA sequence to RNA"""
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
    """Translate RNA sequence into amino acid sequence"""
    seq = transcribe(seq)
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)] # separate codons
    pep = ""
    for c in codons:
        if len(c) < 3:
            break
        pep += codon_match(c)
    
    return pep

def codon_match(codon):
    """Translates single codon to amino acid

    Args:
        codon (List[str,3]): length 3 nucleo sequence

    Returns:
        char: amino acid single-letter code 
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

def input_names(length, reads, i):
    """Constructs protein/read file paths based on benchmark parameters

    Args:
        length (int): average length of each read
        reads (int): total number of reads per file
        i (int): file number

    Returns:
        Tuple[str, str]: protein file path, read file name
    """
    folder = f"{length}-{reads}"
    base = f"{length}-{reads}-{i}"
    return f"{folder}/{base}-protein.txt", f"{folder}/{base}-read.txt"
