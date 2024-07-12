import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import MDAnalysis as mda

# add multi?

# folders with avg files (assumes exact names)
xvg_folders = ['xvg_files_OFF_HMR', 'xvg_files_macrog', 'xvg_files_slipids']

# angle ranges
positive_range = (30, 90)
negative_range = (-90, -30)

# dictionary to store percentages for each folder
percentages_dict = {folder: [] for folder in xvg_folders}

# calculate probability percentages within ranges
def calculate_percentage(xvg_file):

    aux = mda.auxiliary.XVG.XVGReader(xvg_file)
    data = []
    for step in aux:
        step_data = step.data
        if step_data.ndim == 1:
            step_data = step_data.reshape(-1, 2)
        data.extend(step_data)
    
    df = pd.DataFrame(data, columns=['A', 'B'])
    
    # sum of probabilities within ranges
    positive_prob = df[(df['A'] >= positive_range[0]) & (df['A'] <= positive_range[1])]['B'].sum()
    negative_prob = df[(df['A'] >= negative_range[0]) & (df['A'] <= negative_range[1])]['B'].sum()
    
    # calculate total probability to normalize
    total_prob = df['B'].sum()
    
    # sum
    combined_prob = positive_prob + negative_prob
    
    # convert sum to percentage
    combined_percentage = (combined_prob / total_prob) * 100 if total_prob != 0 else 0
    
    return combined_percentage

# loop over each XVG folder and calculate percentages
for folder in xvg_folders:
    xvg_files = [os.path.join(folder, f) for f in sorted(os.listdir(folder)) if f.endswith('.xvg')]
    for xvg_file in xvg_files:
        combined_percent = calculate_percentage(xvg_file)
        percentages_dict[folder].append(combined_percent)

# plot
fig, ax = plt.subplots(figsize=(10, 6))

colors = ['#0333b0', 'green', '#ee7f17']            # adjust as necessary
labels = ['OpenFF', 'MacRog', 'Slipids']

for folder, color, label in zip(xvg_folders, colors, labels):
    percentages = percentages_dict[folder]
    x_labels = np.arange(1, len(percentages) + 1)
    mean_percentages = np.mean(percentages)
    se_percentages = np.std(percentages) / np.sqrt(len(percentages))
    ax.errorbar(x_labels, percentages, yerr=se_percentages, fmt='o', color=color, ecolor=color, capsize=5, label=f'{label} (Mean: {mean_percentages:.2f}%, SE: {se_percentages:.2f}%)')

ax.set_xlabel('Bond #')
ax.set_ylabel('% Isomerization')
ax.set_title('sn-1 Gauche Isomerization')
ax.legend()
ax.set_xticks(np.arange(1, len(percentages_dict[xvg_folders[0]]) + 1))
fig.tight_layout()

plt.savefig('sn1_isomerization_comparison.png', dpi=300)                            # IMAGE SAVE

plt.show()

# output the results into XVG file for each folder
for folder in xvg_folders:
    output_xvg = f'sn1_isomerization_{folder}.xvg'
    with open(output_xvg, 'w') as file:
        file.write('# Probability Percentages within Specified Angle Ranges\n')
        file.write('# Bond#  Combined(30 to 90 and -90 to -30)\n')
        percentages = percentages_dict[folder]
        for i, combined_percent in enumerate(percentages, start=1):
            file.write(f"{i}  {combined_percent:.2f}\n")
