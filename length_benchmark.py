import argparse
import timeit

from less_naive import main as less_naive
from naiveKmerIndex import main as naive

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--length', required=True, type=int, help="read length")
    parser.add_argument('-n', '--number', required=True, type=int, help="number of reads")
    parser.add_argument('-f', '--files', required=True, type=int, help="number of files")
    parser.add_argument('-r', '--repeats', type=int, default=10, help="repeats per files")
    parser.add_argument('-k', '--k', required=True, type=int, help="k value to use")
    parser.add_argument('-m', '--mode', required=True, choices=["naive", "shared", "random"], help="which matching mode to use")
    args = parser.parse_args()

    return args

def input_names(length, reads, i):
    folder = f"{length}-{reads}"
    base = f"{length}-{reads}-{i}"
    return f"{folder}/{base}-protein.txt", f"{folder}/{base}-read.txt"

def benchmark(args, matcher):
    for i in range(args.files):
        protein_file, dna_file = input_names(args.length, args.number, i)
        if args.mode == "random":
            f = lambda: matcher(args.k, protein_file, dna_file, output_file=None, randomize=True, repeats=args.repeats)
        else:
            f = lambda: matcher(args.k, protein_file, dna_file, output_file=None, repeats=args.repeats)
        total_time = timeit.Timer(f).timeit(number=1)

        print(total_time/args.repeats)

if __name__ == "__main__":
    args = parse_args()

    matcher = naive if args.mode == "naive" else less_naive
    print(f"Using {args.mode} matching strategy\n")

    benchmark(args, matcher)
