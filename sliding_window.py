import numpy as np
from matplotlib import pyplot as plt
import argparse
import math
import sys
from alignment import alignment_iterator
import alignment


# Read alignment file
def get_alignments(input_file):
  alignments = []
  for i, line in enumerate(input_file):
    tokens = line.strip().split("\t")
    query_label, query_start, query_seq, target_label, target_start, target_seq = tokens
    query_name = query_label.split(";")[0]
    expected_errors = float(query_label.split(";")[1].split("=")[1])
    alignments.append(alignment.Alignment(int(query_start), query_seq, int(target_start), target_seq, query_name, target_label, expected_errors))
  return alignments


def max_mismatch_sliding_window(alignment, sliding_window_size=50):
  running_mismatch_count = [0] * sliding_window_size
  mismatch_count = 0
  i = 0
  max_total = 0
  for position in alignment_iterator(alignment):
    if position.target_nucleotide != position.query_nucleotide:
      mismatch_count += 1
    running_mismatch_count[i] = mismatch_count
    i = (i + 1) % sliding_window_size
    max_total = max(max_total, mismatch_count - running_mismatch_count[i])
  return max_total


def save_max_sliding_mismatches(alignments, output_file):
  max_mismatch_rates = [max_mismatch_sliding_window(a) for a in alignments]

  plt.cla()
  number_of_bins = max(max_mismatch_rates)
  number_of_bins = round(number_of_bins / math.ceil(number_of_bins/50))
  plt.hist(max_mismatch_rates, number_of_bins, log=True)
  plt.ylim(ymin=0.9)
  plt.xlabel('Max mismatches')
  plt.ylabel('Number of sequences')
  plt.tick_params(which='both', direction='out')
  plt.title('Max number of mismatches in any 50 bases')
  plt.grid(True)
  plt.savefig(output_file)


def print_human_readable(a, value):
  print(">Target_start_{}_len_{}_max_mismatches_{}".format(a.target_start, a.target_length, value))
  print(a.target_seq)
  print(">Query_len_{}_{}".format(a.query_length, a.name))
  print(a.query_seq)


def print_formatted(a):
  print("{};ee={:.4f}\t{}\t{}\t{}\t{}\t{}".format(a.name, a.expected_errors, a.query_start, a.query_seq, a.target_label, a.target_start, a.target_seq))


def print_bad_sequences(alignments, max_mismatch_threshold=5, human_readable=False):
  for a in alignments:
    max_mismatches = max_mismatch_sliding_window(a)
    if max_mismatches > max_mismatch_threshold and human_readable:
      print_human_readable(a, max_mismatches)
    elif max_mismatches > max_mismatch_threshold:
      print_formatted(a)


def print_good_sequences(alignments, max_mismatch_threshold=5, human_readable=False):
  for a in alignments:
    max_mismatches = max_mismatch_sliding_window(a)
    if max_mismatches <= max_mismatch_threshold and human_readable:
      print_human_readable(a, max_mismatches)
    elif max_mismatches <= max_mismatch_threshold:
      print_formatted(a)


def run_from_command_line():
  parser = argparse.ArgumentParser(description="Simple visualisation of fasta/fastq sequences converting bases to color pixels.")
  parser.add_argument('input', type=argparse.FileType('r'), help="the location of the fasta/fastq file to visualise")
  parser.add_argument('-m', '--max-sliding-mismatch-file', type=str, help="the output file for the sliding mismatch visualisation. File format of image is inferred from extension")
  parser.add_argument('-g', '--print-good-sequences', type=int, help="the max number of mismatches in a sliding window of 100 that still counts as good")
  parser.add_argument('-b', '--print-bad-sequences', type=int, help="the max number of mismatches in a sliding window of 100 that still counts as good")
  parser.add_argument('--human-readable', action='store_true', help="human readable")
  args = parser.parse_args()
  alignments = get_alignments(args.input)
  if args.max_sliding_mismatch_file:
    print("Saving sliding mismatch plot " + args.max_sliding_mismatch_file, file=sys.stderr)
    save_max_sliding_mismatches(alignments, args.max_sliding_mismatch_file)
  if args.print_good_sequences:
    print_good_sequences(alignments, args.print_good_sequences, args.human_readable)
  if args.print_bad_sequences:
    print_bad_sequences(alignments, args.print_bad_sequences, args.human_readable)

if __name__ == "__main__":
  run_from_command_line()
