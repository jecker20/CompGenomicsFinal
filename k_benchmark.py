import argparse
import timeit

from common import input_names
from less_naive import main as less_naive
from naiveKmerIndex import main as naive

def parse_args():
    """Parses command-line arguments for the benchmark utility"""
    parser = argparse.ArgumentParser(
        prog="k-mer k-value benchmark",
        description="Run benchmark on a given matching method, over a range of k-values"
    )
    parser.add_argument('-l', '--length', required=True, type=int, help="read length")
    parser.add_argument('-n', '--number', required=True, type=int, help="number of reads")
    parser.add_argument('-f', '--files', required=True, type=int, help="number of files")
    parser.add_argument('-r', '--repeats', type=int, default=10, help="repeats per files")
    parser.add_argument('--kmin', type=int, default=3, help="smallest k value to try")
    parser.add_argument('--kmax', type=int, default=50, help="largest k value to try")
    parser.add_argument('--kstep', type=int, default=2, help="step to use between k values")
    parser.add_argument('-m', '--mode', required=True, choices=["naive", "shared", "random"], help="which matching algorithm to use")
    args = parser.parse_args()

    return args

def benchmark(args, matcher):
    """Run benchmark given parameters

    Args:
        args: command-line arguments
        matcher: protein-read matching function to benchmark

    Returns:
        float: average processing time across all files
    """
    times = []
    for i in range(args.files):
        protein_file, dna_file = input_names(args.length, args.number, i)
        if args.mode == "random":
            f = lambda: matcher(args.k, protein_file, dna_file, output_file=None, randomize=True, repeats=args.repeats)
        else:
            f = lambda: matcher(args.k, protein_file, dna_file, output_file=None, repeats=args.repeats)
        total_time = timeit.Timer(f).timeit(number=1)

        times.append(total_time/args.repeats)

    return sum(times)/args.files

if __name__ == "__main__":
    args = parse_args()

    matcher = naive if args.mode == "naive" else less_naive
    print(f"Using {args.mode} matching strategy\n")

    print("k\ttime")

    for k in range(args.kmin, args.kmax, args.kstep):
        args.k = k
        avg_time = benchmark(args, matcher)
        print(f"{k}\t{avg_time}")
