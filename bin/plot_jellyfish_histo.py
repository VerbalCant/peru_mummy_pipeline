#!/usr/bin/env python
# plot_jellyfish_histo.py
# /u/VerbalCant
# individual step to plot a jellyfish histogram of k-mer frequency
# to a png for display in... wherever you want to display it.

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def plot_jellyfish_histogram(input_file, output_file, k_min=None, bin_size=1):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    data = [line.strip().split() for line in lines]
    df = pd.DataFrame(data, columns=['Count', 'Frequency'])
    df = df.astype({'Count': int, 'Frequency': int})

    if k_min:
        df = df[df['Count'] >= k_min]

    if bin_size != 1:
        df['Binned_Count'] = (df['Count'] // bin_size) * bin_size
        df = df.groupby('Binned_Count').Frequency.sum().reset_index()
        df.rename(columns={'Binned_Count': 'Count'}, inplace=True)

    plt.figure(figsize=(15, 6))
    sns.lineplot(x='Count', y='Frequency', data=df[df['Count'] < 100], marker='o')
    plt.title('Jellyfish K-mer Histogram')
    plt.xlabel('K-mer Count')
    plt.ylabel('Frequency')
    plt.grid(True)

    plt.savefig(output_file)

parser = argparse.ArgumentParser(description='Plot Jellyfish K-mer Histogram.')
parser.add_argument('input_file', type=str, help='Path to the Jellyfish histogram file.')
parser.add_argument('-o', '--output', type=str, default='jellyfish_histogram.png', help='Output PNG file name.')
parser.add_argument('--k-min', type=int, help='Minimum k-mer count to include in the plot.')
parser.add_argument('--bin-size', type=int, default=1, help='Bin size for grouping k-mer counts.')

args = parser.parse_args()
plot_jellyfish_histogram(args.input_file, args.output, args.k_min, args.bin_size)

