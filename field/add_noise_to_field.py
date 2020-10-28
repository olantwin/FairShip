import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from scipy.ndimage import gaussian_filter
import random
import argparse
import os

def plot_my_hist(datum):
    plotData = datum[datum['y'] == 0]
    H, xedges, yedges = np.histogram2d(plotData['x'], plotData['z'], bins=[50, 500], weights=plotData['by'])
    plt.figure(figsize=[20, 10])
    plt.imshow(H, interpolation='nearest', origin='low')
    # plt.colorbar()
    plt.show()

def generate_file(input_fileName, output, xSpace=73, ySpace=128, zSpace=1214, step=2.5, args=None):
    # (min, max, max/stepSize + 1)  in case of Z: (0, nSteps*2.5 - 2.5, nSteps)
    field = pd.read_csv(input_fileName, skiprows=1, sep ='\s+', names=['x', 'y', 'z', 'bx', 'by', 'bz'])

    if args.sidesOnly:
        field_mask[['bx', 'by', 'bz']] = field[['bx', 'by', 'bz']] != 0
        temp_by = np.array(field['by']).reshape([xSpace, ySpace, zSpace])
        temp_by = gaussian_filter(temp_by, sigma=args.sigma)
        field['by'] = temp_by.reshape(-1)
        field['by'] = field['by'] * field_mask['by']
    else:
        field_new = np.zeros(len(field))
        index_range = np.random.randint(0, len(field_new), size=args.nCores)
        field_new[index_range] = random.uniform(-args.peak, args.peak)
        field_new = field_new.reshape([xSpace, ySpace, zSpace])
        field_new = gaussian_filter(field_new, sigma=args.sigma).reshape(-1)
        field_new = field_new / (np.abs(field_new).max()) * (field[['by']] != 0).values.squeeze()
        field['by'] = field['by'] + field['by'] * field_new * args.fraction

    # plot_my_hist(field_mask)
    # plot_my_hist(field)
    # plot_my_hist(field_new)
    # plot_my_hist(rezult)
    field.to_csv(output, sep='\t', header=None, index=None)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process field.')
    parser.add_argument('--input', default=os.path.expandvars("$FAIRSHIP/files/fieldMap.csv"), type=str, action='store')
    parser.add_argument('--output', default=os.path.expandvars("$FAIRSHIP/files/noisy_fieldMap.csv"), type=str, action='store')
    parser.add_argument('--sidesOnly', default=False,  action='store_true')
    parser.add_argument('--sigma', default=30, type=float, action='store')
    parser.add_argument('--nCores', default=1000, type=int, action='store')
    parser.add_argument('--peak', default=500, type=float, action='store')
    parser.add_argument('--fraction', default=0.4, type=float, action='store')
    args = parser.parse_args()

    with open(args.input) as f:
        first_line = f.readline().strip().split()

    generate_file(args.input, args.output,
                  xSpace=int(first_line[0]),
                  ySpace=int(first_line[1]),
                  zSpace=int(first_line[2]),
                  step=float(first_line[3]), args=args)