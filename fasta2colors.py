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

parser = argparse.ArgumentParser(description="Simple visualisation of fasta/fastq sequences converting bases to color pixels.")
parser.add_argument('input', type=argparse.FileType('r'), help="the location of the fasta/fastq file to visualise")
parser.add_argument('output', type=str, help="the output file for the visualisation. File format of image is inferred from extension")
parser.add_argument('--case-sensitive', action='store_true', help="do not ignore case of chars in input file")
parser.add_argument('-s', '--sort', action='store_true', help="sort sequences by sequence length")
parser.add_argument('-c', '--color-palette', type=str, choices=["fixed", "dynamic"], default="fixed")
args = parser.parse_args()

#TODO:
# Blended colors for combo letters

# This color table is sourced from https://github.com/trident01/BioExt-1/blob/master/AlignmentImage.java
fixed_color_table = {
  "A": [255, 0, 0],
  "C": [255, 255, 0],
  "T": [0, 255, 0],
  "G": [190, 0, 95]}

def iterate_fasta(input_file, new_seq_func, new_line_of_seq_func):
  input_file.seek(0)
  ignore_lines = False
  for line in input_file:
    line = line.strip()
    if line[0] == ">": # FASTA Label
      new_seq_func()
    elif line[0] == "+": # FASTQ beginning of quality info
      ignore_lines = True
    elif line[0] == "@": # FASTQ Label
      new_seq_func()
      ignore_lines = False # We are at a new sequence
    elif not ignore_lines:
      new_line_of_seq_func(line)

def fill_color_matrix(input_file, data_matrix, color_map, row_map=None):
  row = -1
  column = 0
  def new_seq_func():
    nonlocal row, column
    row += 1
    column = 0
  def new_line_in_seq_func(line):
    nonlocal row, column, data_matrix, color_map, row_map
    for letter in line:
      colors = color_map(letter)
      for color, value in enumerate(colors):
        if not row_map is None:
          data_matrix[row_map[row], column, color] = value
        else:
          data_matrix[row, column, color] = value
      column += 1
  iterate_fasta(input_file, new_seq_func, new_line_in_seq_func)
 
def get_row_lengths(input_file):
  row_lengths = []
  def new_seq_func():
    nonlocal row_lengths
    row_lengths.append(0)
  def new_line_in_seq_func(line):
    nonlocal row_lengths
    row_lengths[-1] += len(line)
  iterate_fasta(input_file, new_seq_func, new_line_in_seq_func)
  return row_lengths

def file_dimensions(row_lengths):
  return (len(row_lengths), max(row_lengths))

def get_unique_letters(input_file, case_sensitive):
  letters = set()
  def new_line_in_seq_func(line):
    nonlocal letters, case_sensitive
    letters.update(list(line if case_sensitive else line.upper()))
  iterate_fasta(input_file, lambda: None, new_line_in_seq_func)
  letters.discard('\n')
  letters.discard('-')
  return letters 

def make_dynamic_colors(letters):
  import colorsys
  width = 1.0 / len(letters)
  color_palette = {}
  for i, letter in enumerate(letters):
    color_palette[letter] = colorsys.hsv_to_rgb(i * width, 1.0, 255) 
  return color_palette

if args.color_palette == "fixed":
  color_palette = fixed_color_table
elif args.color_palette == "dynamic":
  letters = get_unique_letters(args.input, args.case_sensitive)
  print("Detected dictionary: ", letters)
  color_palette = make_dynamic_colors(letters)

row_lengths = get_row_lengths(args.input)
row_map = None
if args.sort:
  row_map = np.argsort(np.argsort(row_lengths))
dimensions = file_dimensions(row_lengths)
print("Dimensions: " + str(dimensions))
data_matrix = np.full((dimensions[0], dimensions[1], 3), 255, dtype=np.uint8)
fill_color_matrix(args.input, data_matrix, lambda x: color_palette.get(x if args.case_sensitive else x.upper(), [0, 0, 0]), row_map=row_map)
img = scipy.misc.toimage(data_matrix)
img.save(args.output)
