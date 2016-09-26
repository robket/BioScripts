import numpy as np
import scipy.misc
import argparse
from matplotlib import pyplot as plt
import fileinput

parser = argparse.ArgumentParser(description="Make histogram from input data")
parser.add_argument('output', type=str, help="the output file for the histogram. File format of image is inferred from extension")
args = parser.parse_args()

values = []
for line in fileinput.input():
  values.append(float(line.strip()))

plt.cla()
plt.hist(values, 50)
plt.ylim(ymin=0.9)
plt.xlabel('Expected number of errors')
plt.ylabel('Number of sequences')
plt.tick_params(which='both', direction='out')
plt.title('Mismatch Rates')
plt.grid(True)
plt.savefig(output_file)
