import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from io import StringIO

def read_xvg(file):
    with open(file, 'r') as f:
        # skip comments
        lines = [line for line in f if not line.startswith(('#', '@'))]
    # write xvg columns into dataframe
    df = pd.read_csv(StringIO('\n'.join(lines)), delim_whitespace=True, header=None)
    df.columns = ['angle', 'probability']
    return df

def plot_heatmap(correlation_matrix, xticklabels=None, yticklabels=None, xlabel='', ylabel='', title=''):
    plt.figure(figsize=(10, 8))
    mask = np.triu(np.ones_like(correlation_matrix, dtype=bool), k=1)

    sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', vmin=-1, vmax=1, cbar=True,
                xticklabels=xticklabels, yticklabels=yticklabels, mask=mask)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    plt.savefig('TOTAL_CORRELATION_MATRIX.png', dpi=300)        # IMAGE SAVE

    plt.show()

def analyze_dihedrals(prefixes, suffixes, forcefield_names=None, xlabel='Forcefield', ylabel='Forcefield', title=''):
    all_correlations = []

    for prefix in prefixes:
        files = [f"{prefix}{suffix}.xvg" for suffix in suffixes]

        # dict to store xvg data by rows
        data_frames = {}

        # read xvg and store in dict
        for file, forcefield in zip(files, forcefield_names):
            df = read_xvg(file)
            data_frames[forcefield] = df['probability']  # only select second column

        # df with second column
        df = pd.DataFrame(data_frames)

        # pearson correlation coefficient matrix
        correlation_matrix = df.corr(method='pearson')
        all_correlations.append(correlation_matrix)

    # average the correlation matrices to build total
    total_correlation_matrix = sum(all_correlations) / len(all_correlations)

    # plot single matrix
    plot_heatmap(total_correlation_matrix, xticklabels=forcefield_names, yticklabels=forcefield_names,
                 xlabel=xlabel, ylabel=ylabel, title=title)

# torsion ids and main call
forcefield_names = ['OpenFF', 'MacRog', 'Slipids']
prefixes = ['P1_O1_C2_H2', 'P1_O4_C3_H4', 'C2_O1_P1_O4', 'C2_O1_P1_O3', 'C1_C2_O1_P1', 'C2_O1_P1_O2', 'P1_O1_C2_H3', 'P1_O4_C3_C4', 'P1_O4_C3_H5']
suffixes = ['_OFF', '_macrog', '_slipids']

analyze_dihedrals(prefixes, suffixes, forcefield_names, xlabel='Force Field', ylabel='Force Field', title="Total Correlation Matrix of Phosphate Torsions")
