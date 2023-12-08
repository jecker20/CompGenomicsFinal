"""Shared k-mer matching algorithm implementation"""

import argparse
import random

from common import parse_syntheticProt, parse_syntheticDNA, make_orfs, translate

class KMerIndex:
    """k-mer index, capable of handling multiple input sequences"""

    def __init__(self, k, sequences, randomize=False):
        """Constructs index for given sequences and k-value.
        
        Sequences should be a list of tuples, where each tuple is (sequence name, sequence).

        List of all k-mers is created, either sorted by total number of occurences (in decreasing order),
        or in a random order if the `randomize` flag is specified.
        """
        self.k = k
        self.data = sequences
        self.index = {}
        self.counts = {}

        for name, seq in sequences.items():
            self._add_sequence(name, seq)

        self.kmers = list(self.index.keys())
        if randomize:
            random.shuffle(self.kmers)
        else:
            self.kmers.sort(key = lambda s: self.counts[s], reverse=True)

    def _add_sequence(self, name, seq):
        """Internal function to index a new sequence"""
        if self.k == 0 or len(seq) < self.k:
            return
        
        for i in range(len(seq) - self.k + 1):
            s = seq[i:i+self.k]
            if s not in self.index:
                self.index[s] = []
                self.counts[s] = 0
            self.index[s].append((name, i))
            self.counts[s] += 1

    def size(self):
        """Total number of distinct k-mers"""
        return len(self.index)

    def get(self, kmer):
        """Retrieve list of occurences of given k-mer.
        
        Return value is list of tuples in form (sequence name, index of occurence)
        """
        try:
            return self.index[kmer]
        except:
            return []

    def contains(self, kmer):
        """Whether the index contains any occurences of the given k-mer."""
        return kmer in self.index

def is_match(read_hit, protein_hit, reads, proteins):
    """Determines whether (read occurence, protein occurence) pair is a match.
    
    read_hit, protein_hit should be (name, index) format of KMerIndex occurences.
    reads, proteins should be Dict[sequence name, sequence value].
    """
    protein_name, h_loc = protein_hit
    read_name, r_loc = read_hit

    # Given critical assumption that read errors are edit errors only, we can eliminate
    # any hits that *don't* occur at the same location. If we wanted to remove this
    # assumption, we would instead have a tolerance for how far apart the hits can be.
    if h_loc != r_loc:
        return False

    read = reads[read_name]
    protein = proteins[protein_name]

    if read == protein:
        return True

def find_matches(read_index: KMerIndex, protein_index: KMerIndex):
    """Finds all occurences of reads within proteins, given constructed
    k-mer index for both.

    Returns:
        List[Tuple[str, str]]: List of (read name, protein name) matches
    """
    matches = set()
    matched = set()

    for kmer in read_index.kmers:
        p_hits = protein_index.get(kmer)
        r_hits = read_index.get(kmer)
        for p in p_hits:
            for r in r_hits:
                read_name = r[0]
                protein_name = p[0]
                # skips reads we already found a match for
                if read_name in matched:
                    continue

                if is_match(r, p, read_index.data, protein_index.data):
                    matches.add((read_name, protein_name))
                    matched.add(read_name)

    return matches

def main(k, protein_file, read_file, output_file=None, randomize=False, repeats=1):
    """Run our matching algorithm

    Args:
        k (int): k-value for k-mer indices
        protein_file (str): proteins file path (FASTQ-inspired format)
        read_file (str): DNA reads file path (FASTQ-inspired format)
        output_file (str, optional): File path to save matches to. Defaults to None.
        randomize (bool, optional): Whether to randomize k-mer search order. Defaults to False.
        repeats (int, optional): Repeat matching process (for benchmarking). Defaults to 1.
    """
    proteins = parse_syntheticProt(protein_file)
    reads = parse_syntheticDNA(read_file)

    # build the protein index a single time
    protein_index = KMerIndex(k, proteins)

    # re-build read index every trial
    for _ in range(repeats):
        frame_name = lambda read_name, frame_idx: f"{read_name}-{frame_idx+1}"
        peptide_reads = { frame_name(read_name, i) : translate(frame)
                for read_name, read in reads.items()
                for i, frame in enumerate(make_orfs(read)) }

        read_mkers = KMerIndex(k, peptide_reads, randomize)
        matches = find_matches(read_mkers, protein_index)

    if output_file:
        with open(output_file, 'w') as f:
            for read_frame_name, protein_name in matches:
                f.write(f"{read_frame_name},{protein_name}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', '--k', type=int, default=4)
    parser.add_argument('-p', '--proteins', required=True)
    parser.add_argument('-r', '--reads', required=True)
    parser.add_argument('-o', '--output', required=True)    
    args = parser.parse_args()

    main(args.k, args.proteins, args.reads, args.output)
