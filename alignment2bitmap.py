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

parser = argparse.ArgumentParser(description="Simple visualisation of sequence alignments.")
parser.add_argument('input', type=argparse.FileType('r'), help="the location of the alignment (.TSV) file to visualise")
parser.add_argument('output', type=str, help="the output file for the visualisation. File format of image is inferred from extension")
parser.add_argument('-s', '--sort-by', type=str, choices=["mean", "start", "end", "length", "none"], default="mean", help="what parameter is used to sort the sequences")
parser.add_argument('--do-not-crop', action='store_true', help="do not crop the image to the position of the leftmost alignment")
args = parser.parse_args()

def get_start_and_end_points(input_file):
  coords = []
  for line in input_file:
    length, start, end = list(map(int, line.strip().split())) 
    coords.append((start, end))
  return coords

coords = get_start_and_end_points(args.input)
if args.do_not_crop:
  minimum = 0
else:
  minimum = min(coords, key=lambda x: x[0])[0]
maximum = max(coords, key=lambda x: x[1])[1]
dimensions = (len(coords), maximum - minimum)
data_matrix = np.full((dimensions[0], dimensions[1] + 1), 255, dtype=np.uint8)
if args.sort_by == "mean":
  coords.sort(key=sum)
if args.sort_by == "start":
  coords.sort()
elif args.sort_by == "end":
  coords.sort(key=lambda x: (x[1], x[0]))
elif args.sort_by == "length":
  coords.sort(key=lambda x: (x[1] - x[0], x[0]))

for i, (start, end) in enumerate(coords):
  np.put(data_matrix[i], range(start - minimum, end - minimum), 0)
img = scipy.misc.toimage(data_matrix)
img.save(args.output)
