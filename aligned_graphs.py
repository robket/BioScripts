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


def run_from_command_line():
  parser = argparse.ArgumentParser(description="Simple visualisation of fasta/fastq sequences converting bases to color pixels.")
  parser.add_argument('input', type=argparse.FileType('r'), help="the location of the fasta/fastq file to visualise")
  parser.add_argument('-a', '--alignment-output-file', type=str, help="the output file for the alignment map. File format of image is inferred from extension")
  parser.add_argument('-c', '--coverage-output-file', type=str, help="the output file for the coverage visualisation. File format of image is inferred from extension")
  parser.add_argument('-d', '--deletion-dist-file', type=str, help="the output file for the deletion distribution visualisation. File format of image is inferred from extension")
  parser.add_argument('-i', '--insertion-dist-file', type=str, help="the output file for the insertion distribution visualisation. File format of image is inferred from extension")
  parser.add_argument('-m', '--mismatch-rates-file', type=str, help="the output file for the mismatch rate distribution visualisation. File format of image is inferred from extension")
  parser.add_argument('-n', '--nucleotide-output-file', type=str, help="the output file for the nucleotide visualisation. File format of image is inferred from extension")
  parser.add_argument('-M', '--n-mismatched', type=int, help="print the n most mismatched sequences")
  parser.add_argument('--case-sensitive', action='store_true', help="do not ignore case of chars in input file")
  args = parser.parse_args()

  alignments = get_alignments(args.input)
  if args.insertion_dist_file:
    print("Saving insertion distribution plot " + args.insertion_dist_file, file=sys.stderr)
    alignment.save_insertion_or_deletion_dist(alignments, args.insertion_dist_file)
  if args.deletion_dist_file:
    print("Saving deletion distribution plot " + args.deletion_dist_file, file=sys.stderr)
    alignment.save_insertion_or_deletion_dist(alignments, args.deletion_dist_file, insertion_not_deletion=False)
  if args.mismatch_rates_file:
    print("Saving mismatch rates plot " + args.mismatch_rates_file, file=sys.stderr)
    alignment.save_mismatch_rates(alignments, args.mismatch_rates_file, ignore_case=(not args.case_sensitive))
  if args.coverage_output_file:
    print("Saving coverage map to " + args.coverage_output_file, file=sys.stderr)
    alignment.save_coverage_map(alignments, args.coverage_output_file)
  if args.nucleotide_output_file:
    print("Saving nucleotide map to " + args.nucleotide_output_file, file=sys.stderr)
    alignment.save_nucleotide_map(alignments, args.nucleotide_output_file, ignore_case=(not args.case_sensitive))
  if args.n_mismatched:
    print("{} most mismatched sequences:".format(args.n_mismatched), file=sys.stderr)
    alignment.print_most_mismatched(alignments, args.n_mismatched)
  if args.alignment_output_file:
    print("Saving alignment map to " + args.alignment_output_file, file=sys.stderr)
    coords = [(a.target_start, a.target_start + a.target_length) for a in alignments]
    alignment.save_alignment_map(coords, args.alignment_output_file)

if __name__ == "__main__":
  run_from_command_line()
