import argparse
import sys

from kmer_index import KMerIndex

from common import parse_syntheticProt, parse_syntheticDNA, make_orfs, translate

def filter_out(hit, protein_index):
    name, location = hit
    length = len(protein_index.data)
    return False # TODO: can we filter out for all reads simultaneously?

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

def main(k, protein_file, read_file, output_file):
    synProts = parse_syntheticProt(protein_file)
    reads = parse_syntheticDNA(read_file)

    frame_name = lambda read_name, frame_idx: f"{read_name}-{frame_idx+1}"
    peptide_reads = { frame_name(read_name, i) : translate(frame)
               for read_name, read in reads.items()
               for i, frame in enumerate(make_orfs(read)) }
    
    proteins = synProts

    read_mkers = KMerIndex(k, peptide_reads)
    protein_index = KMerIndex(k, proteins)

    matches = find_matches(read_mkers, protein_index)

    #with open(output_file, 'w') as f:
    #    for read_frame_name, protein_name in matches:
    #        f.write(f"{read_frame_name},{protein_name}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', '--k', type=int, default=4)
    parser.add_argument('-p', '--proteins', required=True)
    parser.add_argument('-r', '--reads', required=True)
    parser.add_argument('-o', '--output', required=True)    
    args = parser.parse_args()

    main(args.k, args.proteins, args.reads, args.output)
