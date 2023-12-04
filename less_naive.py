import argparse

from common import parse_syntheticProt, parse_syntheticDNA, make_orfs, translate

class KMerIndex:
    def __init__(self, k, sequences):
        self.k = k
        self.data = sequences
        self.index = {}
        self.counts = {}

        for name, seq in sequences.items():
            self.add_sequence(name, seq)

        self.most_shared = list(self.index.keys())
        self.most_shared.sort(key = lambda s: self.counts[s], reverse=True)

    def add_sequence(self, name, seq):
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
        return len(self.index)

    def get(self, key):
        try:
            return self.index[key]
        except:
            return []

    def contains(self, key):
        return key in self.index

def is_match(r_hit, p_hit, k, reads, proteins):
    protein_name, h_loc = p_hit
    read_name, r_loc = r_hit

    if h_loc != r_loc:
        return False

    read = reads[read_name]
    protein = proteins[protein_name]

    if read == protein:
        return True

def find_matches(read_kmers: KMerIndex, protein_kmers: KMerIndex):
    matches = set()
    matched = set()

    for kmer in read_kmers.most_shared:
        p_hits = protein_kmers.get(kmer)
        r_hits = read_kmers.get(kmer)
        for p in p_hits:
            for r in r_hits:
                read_name = r[0]
                protein_name = p[0]
                # skips reads we already found a match for
                if read_name in matched:
                    continue

                if is_match(r, p, read_kmers.k, read_kmers.data, protein_kmers.data):
                    matches.add((read_name, protein_name))
                    matched.add(read_name)

    return matches

def main(k, protein_file, read_file, output_file=None, repeats=1):
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

        read_mkers = KMerIndex(k, peptide_reads)
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
