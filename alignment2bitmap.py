try:
  import numpy as np
  import scipy.misc
except ImportError:
 print("Scipy does not seem to be available. Try 'pip install scipy'\nError:")
 raise
try:
  from scipy.misc import toimage
except ImportError:
  print("Python Image Library does not seem to be available. Try 'pip install image'\nError:")
  raise
import argparse
import alignment

def get_start_and_end_points(input_file):
  coords = []
  for line in input_file:
    length, start, end = list(map(int, line.strip().split())) 
    coords.append((start, end))
  return coords


def run_from_command_line():
  parser = argparse.ArgumentParser(description="Simple visualisation of sequence alignments.")
  parser.add_argument('input', type=argparse.FileType('r'), help="the location of the alignment (.TSV) file to visualise")
  parser.add_argument('output', type=str, help="the output file for the visualisation. File format of image is inferred from extension")
  parser.add_argument('-s', '--sort-by', type=str, choices=["mean", "start", "end", "length", "none"], default="mean", help="what parameter is used to sort the sequences")
  parser.add_argument('--do-not-crop', action='store_true', help="do not crop the image to the position of the leftmost alignment")
  parser.add_argument('--no-text', action='store_true', help="do display initial offset in the image")
  args = parser.parse_args()

  coords = get_start_and_end_points(args.input)
  sort_key = None
  if args.sort_by == "mean":
    sort_key = sum
  if args.sort_by == "start":
    sort_key = lambda x: x
  elif args.sort_by == "end":
    sort_key = lambda x: (x[1], x[0])
  elif args.sort_by == "length":
    sort_key = lambda x: (x[1] - x[0], x[0])
  alignment.save_alignment_map(coords, args.output, sort_key=sort_key, crop=not args.do_not_crop, no_text=args.no_text)

if __name__ == "__main__":
  run_from_command_line()
