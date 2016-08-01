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
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw


def get_start_and_end_points(input_file):
  coords = []
  for line in input_file:
    length, start, end = list(map(int, line.strip().split())) 
    coords.append((start, end))
  return coords


def save_alignments(coords, output_file, sort_key=sum, crop=True, no_text=False):
  if crop:
    minimum = min(coords, key=lambda x: x[0])[0]
  else:
    minimum = 0
  maximum = max(coords, key=lambda x: x[1])[1]
  dimensions = (len(coords), maximum - minimum)
  data_matrix = np.full((dimensions[0], dimensions[1] + 1), 255, dtype=np.uint8)
  if sort_key is not None:
    coords.sort(key=sort_key)
  for i, (start, end) in enumerate(coords):
    # np.put(data_matrix[i], range(start - minimum, end - minimum), 0)
    data_matrix[i, (start - minimum):(end - minimum)] = 0
  img = scipy.misc.toimage(data_matrix)
  if minimum > 0 and not no_text:
    draw = ImageDraw.Draw(img)
    font = ImageFont.truetype("OpenSans-Regular.ttf", int(dimensions[0] / 40))
    text = "Offset: " + str(minimum)
    draw.text((int(dimensions[1] * 0.95 - font.getsize(text)[0]), int(dimensions[1] * 0.05)), text, font=font)
  img.save(output_file)


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
  save_alignments(coords, args.output, sort_key=sort_key, crop=not args.do_not_crop, no_text=args.no_text)

if __name__ == "__main__":
  run_from_command_line()
