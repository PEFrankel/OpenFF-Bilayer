import os
import numpy as np
import json
import pandas as pd
import MDAnalysis as mda

"""
This script analyzes dihedral distribution XVGs to identify rotamer conformations and 
sequences in lipid acyl chains according to these definitions:
- trans (t): φ < -150° or φ > 150°
- gauche- (g-): -90° ≤ φ < -30°
- gauche+ (g+): 30° < φ ≤ 90°

The script calculates:
- eg (end gauche): gauche at terminal positions
- gg (double gauche): consecutive gauche conformations
- gtg: g+tg+ or g-tg- sequences
- gtg' (kinks): g+tg- or g-tg+ sequences
- ng: total number of gauche conformations

Dihedral distribution XVGs of the sn-1 chain for POPC lipids (SMILES string specific for OpenFF) 
can be obtained from "dihedral_distributions.py"
"""

xvg_folders = ['xvg_files_OFF', 'xvg_files_charmm36', 'xvg_files_lipid21']

trans_range = [(-180, -150), (150, 180)]
gauche_minus_range = (-90, -30)
gauche_plus_range = (30, 90)

results = {folder: {} for folder in xvg_folders}

def identify_rotamer(angle):
    if (angle < -150) or (angle > 150): #t
        return 't'
    elif -90 <= angle < -30: #g-
        return 'g-'
    elif 30 < angle <= 90: #g+
        return 'g+'
    else:
        return 'other'

def get_angle_data(xvg_file):
    aux = mda.auxiliary.XVG.XVGReader(xvg_file)
    data = []
    for step in aux:
        step_data = step.data
        if step_data.ndim == 1:
            step_data = step_data.reshape(-1, 2)
        data.extend(step_data)
    
    df = pd.DataFrame(data, columns=['Angle', 'Probability'])
    return df

def calculate_rotamer_distributions(xvg_files):
    avg_trans = []
    avg_gauche_minus = []
    avg_gauche_plus = []
    
    for xvg_file in xvg_files:
        df = get_angle_data(xvg_file)
        
        trans_prob = 0
        for t_range in trans_range:
            trans_prob += df[(df['Angle'] >= t_range[0]) & (df['Angle'] <= t_range[1])]['Probability'].sum()
        
        gauche_minus_prob = df[(df['Angle'] >= gauche_minus_range[0]) & 
                              (df['Angle'] <= gauche_minus_range[1])]['Probability'].sum()
        
        gauche_plus_prob = df[(df['Angle'] >= gauche_plus_range[0]) & 
                             (df['Angle'] <= gauche_plus_range[1])]['Probability'].sum()
        
        total_prob = df['Probability'].sum()
        
        trans_prob = trans_prob / total_prob if total_prob != 0 else 0
        gauche_minus_prob = gauche_minus_prob / total_prob if total_prob != 0 else 0
        gauche_plus_prob = gauche_plus_prob / total_prob if total_prob != 0 else 0
        
        avg_trans.append(trans_prob)
        avg_gauche_minus.append(gauche_minus_prob)
        avg_gauche_plus.append(gauche_plus_prob)
    
    return avg_trans, avg_gauche_minus, avg_gauche_plus

def calculate_rotamer_sequences(avg_trans, avg_gauche_minus, avg_gauche_plus):
    n_positions = len(avg_trans)
    gauche_any = []
    
    for i in range(n_positions):
        gauche_any.append(avg_gauche_minus[i] + avg_gauche_plus[i])
    
    eg = 0 #first and last positions
    if n_positions > 0:
        eg = (gauche_any[0] + gauche_any[-1]) / 2
    
    gg = 0 #consecutive gauche conformations
    for i in range(n_positions - 1):
        gg += gauche_any[i] * gauche_any[i+1]
    
    gtg = 0 # g+tg+ or g-tg-
    gtg_prime = 0 # g+tg- or g-tg+ (kinks)
    
    for i in range(n_positions - 2):
        gtg += avg_gauche_plus[i] * avg_trans[i+1] * avg_gauche_plus[i+2]
        gtg += avg_gauche_minus[i] * avg_trans[i+1] * avg_gauche_minus[i+2]
        
        #kinks
        gtg_prime += avg_gauche_plus[i] * avg_trans[i+1] * avg_gauche_minus[i+2]
        gtg_prime += avg_gauche_minus[i] * avg_trans[i+1] * avg_gauche_plus[i+2]
    
    gtg_total = gtg + gtg_prime
    
    ng = sum(gauche_any)
    
    return {
        'eg': eg,
        'gg': gg,
        'gtg': gtg,
        'gtg_prime': gtg_prime,
        'gtg_total': gtg_total,
        'ng': ng
    }

for folder in xvg_folders:
    xvg_files = [os.path.join(folder, f) for f in sorted(os.listdir(folder)) if f.endswith('.xvg')]
    avg_trans, avg_gauche_minus, avg_gauche_plus = calculate_rotamer_distributions(xvg_files)
    sequences = calculate_rotamer_sequences(avg_trans, avg_gauche_minus, avg_gauche_plus)
    results[folder] = sequences

formatted_results = {}
for folder, folder_results in results.items():
    system_name = folder.replace('xvg_files_', '')
    formatted_results[system_name] = {
        'eg': folder_results['eg'],
        'gg': folder_results['gg'],
        'gtg_prime': folder_results['gtg_prime'],
        'gtg_total': folder_results['gtg_total'],
        'ng': folder_results['ng']
    }

with open('rotamer_analysis_results.json', 'w') as json_file:
    json.dump(formatted_results, json_file, indent=4)