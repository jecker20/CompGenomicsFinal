import argparse
import timeit

from common import input_names
from less_naive import main as less_naive
from naiveKmerIndex import main as naive

def parse_args():
    """Parses command-line arguments for the benchmark utility"""

    parser = argparse.ArgumentParser(
        prog="matching benchmark",
        description="Run benchmark on a given matching method"
    )
    parser.add_argument('-l', '--length', required=True, type=int, help="read length")
    parser.add_argument('-n', '--number', required=True, type=int, help="number of reads")
    parser.add_argument('-f', '--files', required=True, type=int, help="number of files")
    parser.add_argument('-r', '--repeats', type=int, default=10, help="repeats per files")
    parser.add_argument('-k', '--k', required=True, type=int, help="k value to use")
    parser.add_argument('-m', '--mode', required=True, choices=["naive", "shared", "random"], help="which matching algorithm to use")
    args = parser.parse_args()

    return args

def benchmark(args, matcher):
    """Run benchmark given parameters

    Runnings matching algorith args.repeats times per file,
    and prints the average processing time for each file.

    Args:
        args: command-line arguments
        matcher: protein-read matching function to benchmark
    """
    for i in range(args.files):
        protein_file, dna_file = input_names(args.length, args.number, i)
        if args.mode == "random":
            f = lambda: matcher(args.k, protein_file, dna_file, output_file=None, randomize=True, repeats=args.repeats)
        else:
            f = lambda: matcher(args.k, protein_file, dna_file, output_file=None, repeats=args.repeats)
        total_time = timeit.Timer(f).timeit(number=1)

        print(total_time)

if __name__ == "__main__":
    args = parse_args()

    matcher = naive if args.mode == "naive" else less_naive
    print(f"Using {args.mode} matching strategy\n")

    benchmark(args, matcher)
