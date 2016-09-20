import numpy as np
from matplotlib import pyplot as plt
from scipy.misc import toimage
from collections import defaultdict, Counter
from types import SimpleNamespace
from PIL import ImageFont
from PIL import ImageDraw


# This color table is sourced from https://github.com/trident01/BioExt-1/blob/master/AlignmentImage.java
FIXED_COLOR_TABLE = defaultdict(lambda: [0, 0, 0], {
  "A": [255, 0, 0],
  "C": [255, 255, 0],
  "T": [0, 255, 0],
  "G": [190, 0, 95],
  "-": [150, 150, 150]})


class Alignment:
  def __init__(self, query_start, query_seq, target_start, target_seq):
    self.query_start = int(query_start) - 1
    self.query_seq = query_seq
    query_gap_count = query_seq.count("-")
    self.query_length = len(query_seq) - query_gap_count
    self.target_start = int(target_start) - 1
    self.target_seq = target_seq
    target_gap_count = target_seq.count("-")
    self.target_length = len(target_seq) - target_gap_count
    self.no_gap_length = len(target_seq) - target_gap_count - query_gap_count
    if len(target_seq) != len(query_seq):
      raise ValueError("Length of target sequence not equal to length of query sequence")


def _alignment_iterator(alignment, ignore_case=True, include_gaps=False):
  target_index = 0
  target_offset = 0
  query_index = 0
  while target_index < len(alignment.target_seq) and query_index < len(alignment.query_seq):
    if alignment.target_seq[target_index] == "-": # If it is an insertion
      target_offset += 1
    elif alignment.query_seq[query_index] != "-" or include_gaps:
      reference_index = alignment.target_start + target_index - target_offset
      query_nucleotide = alignment.query_seq[query_index].upper() if ignore_case else alignment.query_seq[query_index]
      target_nucleotide = alignment.target_seq[target_index].upper() if ignore_case else alignment.target_seq[target_index]
      yield SimpleNamespace(reference_index=reference_index,
                            target_nucleotide=target_nucleotide,
                            query_nucleotide=query_nucleotide)
    target_index += 1
    query_index += 1


def count_mismatches(alignment, ignore_case=True):
  mismatch_count = 0
  for position in _alignment_iterator(alignment, ignore_case):
    if position.target_nucleotide != position.query_nucleotide:
      mismatch_count += 1
  return mismatch_count


def save_mismatch_rates(alignments, output_file, ignore_case=True):
  mismatch_rates = [count_mismatches(a, ignore_case) / a.no_gap_length for a in alignments]

  plt.cla()
  plt.hist(mismatch_rates, 50)
  plt.xlabel('Rate of mismatches')
  plt.ylabel('Number of sequences')
  plt.title('Mismatch Rates')
  plt.grid(True)
  plt.savefig(output_file)


def gap_distribution(sequence):
  dist = Counter()
  count_length = 0
  for char in sequence:
    if char == "-":
      count_length += 1
    elif count_length > 0:
      dist[count_length] += 1
      count_length = 0
  if count_length > 0:
    dist[count_length] += 1
  return dist


def save_insertion_or_deletion_dist(alignments, output_file, insertion_not_deletion=True):
  size_counter = Counter()
  for a in alignments:
    size_counter += gap_distribution(a.target_seq if insertion_not_deletion else a.query_seq)

  sizes, counts = zip(*size_counter.items())
  plt.cla()
  n, bins, patches = plt.hist(sizes, 50, weights=counts)
  plt.xlabel('Size of insertion' if insertion_not_deletion else 'Size of deletion')
  plt.ylabel('Count')
  plt.title('Insertion size distribution' if insertion_not_deletion else 'Gap size distribution')
  plt.grid(True)
  plt.savefig(output_file)


# Get nucleotide distribution
def nucleotide_distribution(alignments, ignore_case=False, include_gaps=True):
  max_index = 0
  distribution = defaultdict(Counter)
  for a in alignments:
    for position in _alignment_iterator(a, ignore_case, include_gaps):
      distribution[position.reference_index][position.query_nucleotide] += 1
    max_index = max(max_index, a.target_start + a.target_length)
  return [distribution[i] for i in range(max_index)]


def save_nucleotide_map(alignments, output, ignore_case=True, include_gaps=True):
  nucleotides = nucleotide_distribution(alignments, ignore_case, include_gaps)
  width = len(nucleotides)
  keys = set()
  for distribution_at_base in nucleotides:
    keys.update(set(distribution_at_base.keys()))
  keys = sorted(list(keys), key=lambda x: "ZZZ" if x == "-" else x)
  nucleotide_count_array = np.zeros((len(keys), width), dtype=np.uint32)
  for i, key in enumerate(keys):
    for j, counts in enumerate(nucleotides):
      nucleotide_count_array[i, j] = counts[key]
  cum_sum = nucleotide_count_array.cumsum(axis=0)
  height = cum_sum[-1,].max()
  data_matrix = np.full((height, width, 3), 255, dtype=np.uint8)
  for x in range(width):
    for i, key in enumerate(keys):
      start = 0 if i == 0 else cum_sum[i - 1, x]
      end = cum_sum[i, x]
      data_matrix[start:end, x, 0:3] = FIXED_COLOR_TABLE[key]
  img = toimage(data_matrix[::-1,])
  img.save(output)


# Get coverage map
def coverage_map(alignments, include_gaps=False):
  max_index = 0
  coverage = Counter()
  for a in alignments:
    for position in _alignment_iterator(a, True, include_gaps):
      coverage[position.reference_index] += 1
    max_index = max(max_index, a.target_start + a.target_length)
  return [coverage[i] for i in range(max_index)]


def save_coverage_map(alignments, output):
  coverage_with_gaps = coverage_map(alignments, True)
  coverage_without_gaps = coverage_map(alignments, False)
  width = len(coverage_with_gaps)
  height = max(coverage_with_gaps)
  data_matrix = np.full((height, width, 3), 255, dtype=np.uint8)
  for x in range(width):
    y1 = coverage_without_gaps[x]
    y2 = coverage_with_gaps[x]
    #data_matrix[height-y1:height, x, 0:3] = 0
    #data_matrix[height-y2:height-y1, x, 0:3] = 127
    data_matrix[0:y1, x, 0:3] = 0
    data_matrix[y1:y2, x, 0:3] = 127
  img = toimage(data_matrix[::-1,])
  img.save(output)


def save_alignment_map(coords, output_file, sort_key=sum, crop=True, no_text=False):
  if crop:
    minimum = min(coords, key=lambda x: x[0])[0]
  else:
    minimum = 0
  maximum = max(coords, key=lambda x: x[1])[1]
  dimensions = (len(coords), maximum - minimum)
  data_matrix = np.full((dimensions[0], dimensions[1] + 1), 255, dtype=np.uint8)
  if sort_key is not None:
    coords.sort(key=sort_key)

  is_multiple_alignment = len(coords[0]) > 3 and type(coords[0][3]) == list

  # Greyscale over the bounds (or black if not multiple alignment)
  for i, coord in enumerate(coords):
    start = coord[0]
    end = coord[1]
    # np.put(data_matrix[i], range(start - minimum, end - minimum), 0)
    data_matrix[i, (start - minimum):(end - minimum)] = 196 if is_multiple_alignment else 0

  # Black over the subalignments, if any
  if is_multiple_alignment:
    for i, coord in enumerate(coords):
      for subalignment in coord[3]:
        start = subalignment[0]
        end = subalignment[1]
        # np.put(data_matrix[i], range(start - minimum, end - minimum), 0)
        data_matrix[i, (start - minimum):(end - minimum)] = 0

  img = toimage(data_matrix)

  # Text in the bottom right corner
  if not no_text:
    draw = ImageDraw.Draw(img)
    font = ImageFont.truetype("OpenSans-Regular.ttf", int(dimensions[0] / 40))
    text = "Max: " + str(maximum)
    draw.text((int(dimensions[1] * 0.95 - font.getsize(text)[0]), int(dimensions[0] * 0.09 - font.getsize(text)[1])), text, fill="gray", font=font)
  if minimum > 0 and not no_text:
    draw = ImageDraw.Draw(img)
    font = ImageFont.truetype("OpenSans-Regular.ttf", int(dimensions[0] / 40))
    text = "Offset: " + str(minimum)
    draw.text((int(dimensions[1] * 0.95 - font.getsize(text)[0]), int(dimensions[0] * 0.11)), text, fill="gray", font=font)
  img.save(output_file)
