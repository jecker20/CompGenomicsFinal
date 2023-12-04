import argparse
import timeit

from less_naive import main as less_naive

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--length', required=True, type=int, help="read length")
    parser.add_argument('-n', '--number', required=True, type=int, help="number of reads")
    parser.add_argument('-f', '--files', required=True, type=int, help="number of files")
    parser.add_argument('-r', '--repeats', required=True, type=int, help="repeats per files")
    parser.add_argument('-k', '--k', required=True, type=int, help="k value to use")
    args = parser.parse_args()

    return args

def input_names(length, reads, i):
    base = f"{length}-{reads}/{length}-{reads}-{i}"
    return f"{base}-protein.txt", f"{base}-read.txt"

def benchmark(args):
    times = []
    for i in range(args.files):
        protein_file, dna_file = input_names(args.length, args.number, i)
        f = lambda: less_naive(args.k, protein_file, dna_file, "output_matches.txt")
        total_time = timeit.Timer(f).timeit(number=args.repeats)
        times.append(total_time/args.repeats)

    print(sum(times)/args.files)

if __name__ == "__main__":
    args = parse_args()
    benchmark(args)
