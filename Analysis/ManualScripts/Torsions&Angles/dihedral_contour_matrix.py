import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from io import StringIO

def read_xvg(file):
    with open(file, 'r') as f:
        lines = [line for line in f if not line.startswith(('#', '@'))]
    df = pd.read_csv(StringIO('\n'.join(lines)), delim_whitespace=True, header=None)
    df.columns = ['angle', 'probability']
    return df

def create_custom_matrix(torsions, forcefield, suffix, save_filename=None):
    data_frames = {}
    
    for torsion in torsions:
        file = f"{torsion}{suffix}.xvg"
        df = read_xvg(file)
        data_frames[torsion] = df['probability']
    
    df = pd.DataFrame(data_frames)
    
    correlation_matrix = df.corr(method='pearson')
    
    fig, axes = plt.subplots(len(torsions), len(torsions), figsize=(8, 8), constrained_layout=True)
    
    for i, torsion_i in enumerate(torsions):
        for j, torsion_j in enumerate(torsions):
            ax = axes[i, j]
            
            if i > j:
                # Scatter plot for 'bottom-left' triangle
                ax.scatter(df[torsion_j], df[torsion_i], alpha=0.6, s=4, color='#0333b0')  # Scatter plot color
                if j == 0:
                    ax.set_ylabel(torsion_i, fontsize=6)
                if i == len(torsions) - 1:
                    ax.set_xlabel(torsion_j, fontsize=6)
            elif i == j:
                # Probability distribution plot on 'diagonal'
                angles = np.linspace(-180, 180, len(df[torsion_i]))
                ax.plot(angles, df[torsion_i], linewidth=0.4, color='#0333b0')  # Line plot color
                ax.set_xticks([-180, 0, 180])
                ax.set_xlim([-180, 180])
                ax.set_title(torsion_i, fontsize=8, pad=10)
            else:
                # Correlation value on 'upper-right' triangle
                correlation_value = correlation_matrix.loc[torsion_i, torsion_j]
                
                if correlation_value == 0:
                    color = 'lightgrey'
                else:
                    color = plt.cm.bwr((correlation_value + 1) / 2)  # gradient for other values
                
                ax.patch.set_facecolor(color)
                ax.text(0.5, 0.5, f"{correlation_value:.2f}",
                        horizontalalignment='center', verticalalignment='center',
                        fontsize=8 + abs(correlation_value) * 8,
                        color='black')
                ax.set_xticks([])
                ax.set_yticks([])
            ax.tick_params(axis='both', which='major', labelsize=4)
            ax.tick_params(axis='both', which='minor', labelsize=4)
    
    # Colorbar
    cax = fig.add_axes([0.95, 0.2, 0.02, 0.6])  #left, bottom, width, height
    cmap = plt.cm.bwr
    norm = plt.Normalize(-1, 1)
    cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_label('Correlation', fontsize=10)
    cb.ax.tick_params(labelsize=8)
    
    plt.subplots_adjust(wspace=0.0, hspace=0, right=0.7)
    
    plt.suptitle(f"{forcefield} Force Field", fontsize=12)

    if save_filename:
        plt.savefig(save_filename, dpi=300)
    
    plt.show()

# torsions used
torsions = ['P1_O4_C3_C4','P1_O4_C3_H5','P1_O4_C3_H4','C2_O1_P1_O4','C2_O1_P1_O3','C2_O1_P1_O2','P1_O1_C2_H2','P1_O1_C2_H3','C1_C2_O1_P1']
forcefield = 'Open'
suffix = '_OFF'
save_filename = 'custom_matrix.png'

create_custom_matrix(torsions, forcefield, suffix, save_filename)
