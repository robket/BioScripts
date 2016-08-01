# -userfields qilo+qrow+tilo+trow

try:
  import numpy as np
  import scipy.misc
  from matplotlib import pyplot as plt
except ImportError:
  print("Scipy does not seem to be available. Try 'pip install scipy'\nError:")
  raise
try:
  from scipy.misc import toimage
except ImportError:
  print("Python Image Library does not seem to be available. Try 'pip install image'\nError:")
  raise
import argparse
import sys
import alignment


# Read alignment file
def get_alignments(input_file):
  alignments = []
  for line in input_file:
    query_start, query_seq, target_start, target_seq = line.strip().split("\t")
    alignments.append(alignment.Alignment(int(query_start), query_seq, int(target_start), target_seq))
  return alignments


def print_most(alignments, values, n, keyword="value", percentage=True):
  values = np.array(values)
  indexes = sorted(np.argpartition(-values, n)[:n], key=values.__getitem__, reverse=True)
  for i in indexes:
    a = alignments[i]
    print(">Target_start_{}_len_{}_{}_{:.2f}".format(a.target_start, a.target_length, keyword, values[i] * (100 if percentage else 1)))
    print(a.target_seq)
    print(">Query_len_{}".format(a.query_length))
    print(a.query_seq)


def run_from_command_line():
  parser = argparse.ArgumentParser(description="Simple visualisation of fasta/fastq sequences converting bases to color pixels.")
  parser.add_argument('input', type=argparse.FileType('r'), help="the location of the fasta/fastq file to visualise")
  parser.add_argument('-a', '--alignment-output-file', type=str, help="the output file for the alignment map. File format of image is inferred from extension")
  parser.add_argument('-c', '--coverage-output-file', type=str, help="the output file for the coverage visualisation. File format of image is inferred from extension")
  parser.add_argument('-d', '--deletion-dist-file', type=str, help="the output file for the deletion distribution visualisation. File format of image is inferred from extension")
  parser.add_argument('-i', '--insertion-dist-file', type=str, help="the output file for the insertion distribution visualisation. File format of image is inferred from extension")
  parser.add_argument('-m', '--mismatch-rates-file', type=str, help="the output file for the mismatch rate distribution visualisation. File format of image is inferred from extension")
  parser.add_argument('-n', '--nucleotide-output-file', type=str, help="the output file for the nucleotide visualisation. File format of image is inferred from extension")
  parser.add_argument('-D', '--n-deleted', type=int, help="print the n sequences with the most deletions")
  parser.add_argument('-I', '--n-inserted', type=int, help="print the n sequences with the most insertions")
  parser.add_argument('-M', '--n-mismatched', type=int, help="print the n most mismatched sequences")
  parser.add_argument('--case-sensitive', action='store_true', help="do not ignore case of chars in input file")
  args = parser.parse_args()

  alignments = get_alignments(args.input)
  if args.alignment_output_file:
    print("Saving alignment map to " + args.alignment_output_file, file=sys.stderr)
    coords = [(a.target_start, a.target_start + a.target_length) for a in alignments]
    alignment.save_alignment_map(coords, args.alignment_output_file)
  if args.coverage_output_file:
    print("Saving coverage map to " + args.coverage_output_file, file=sys.stderr)
    alignment.save_coverage_map(alignments, args.coverage_output_file)
  if args.deletion_dist_file:
    print("Saving deletion distribution plot " + args.deletion_dist_file, file=sys.stderr)
    alignment.save_insertion_or_deletion_dist(alignments, args.deletion_dist_file, insertion_not_deletion=False)
  if args.insertion_dist_file:
    print("Saving insertion distribution plot " + args.insertion_dist_file, file=sys.stderr)
    alignment.save_insertion_or_deletion_dist(alignments, args.insertion_dist_file)
  if args.mismatch_rates_file:
    print("Saving mismatch rates plot " + args.mismatch_rates_file, file=sys.stderr)
    alignment.save_mismatch_rates(alignments, args.mismatch_rates_file, ignore_case=(not args.case_sensitive))
  if args.nucleotide_output_file:
    print("Saving nucleotide map to " + args.nucleotide_output_file, file=sys.stderr)
    alignment.save_nucleotide_map(alignments, args.nucleotide_output_file, ignore_case=(not args.case_sensitive))
  if args.n_deleted:
    print("{} sequences with the most pct deletions:".format(args.n_deleted), file=sys.stderr)
    values = [(len(a.query_seq) - a.query_length) / a.target_length for a in alignments]
    print_most(alignments, values, args.n_deleted, "deletionpct")
  if args.n_inserted:
    print("{} sequences with the most insertions:".format(args.n_inserted), file=sys.stderr)
    values = [(len(a.target_seq) - a.target_length) / a.target_length for a in alignments]
    print_most(alignments, values, args.n_inserted, "insertionpct")
  if args.n_mismatched:
    print("{} most mismatched sequences:".format(args.n_mismatched), file=sys.stderr)
    values = [alignment.count_mismatches(a) / a.no_gap_length for a in alignments]
    print_most(alignments, values, args.n_mismatched, "mismatchpct")

if __name__ == "__main__":
  run_from_command_line()
