import json
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

def load_json_data(json_file):
    with open(json_file, 'r') as f:
        return json.load(f)

def parse_data(json_data):
    head_group = {
        'γ': [], 'β': [], 'α': [], 'g1': [], 'g2': [], 'g3': [],
    }
    sn1 = defaultdict(list)
    sn2 = defaultdict(list)

###mapping file keys for chains###

    sn1_keys = [
        'M_G1C3_M', 'M_G1C4_M', 'M_G1C5_M', 'M_G1C6_M', 'M_G1C7_M', 'M_G1C8_M', 'M_G1C9_M',
        'M_G1C10_M', 'M_G1C11_M', 'M_G1C12_M', 'M_G1C13_M', 'M_G1C14_M', 'M_G1C15_M', 'M_G1C16_M', 'M_G1C17_M'
    ]
    
    sn2_keys = [
        'M_G2C3_M', 'M_G2C4_M', 'M_G2C5_M', 'M_G2C6_M', 'M_G2C7_M', 'M_G2C8_M', 'M_G2C9_M', 
        'M_G2C10_M', 'M_G2C11_M', 'M_G2C12_M', 'M_G2C13_M', 'M_G2C14_M', 'M_G2C15_M', 'M_G2C16_M', 'M_G2C17_M', 'M_G2C18_M', 'M_G2C19_M'
    ]
    
    for key, values in json_data.items():
        atom_name = key.split()[0]
        
        if len(values[0]) == 1:
            order_param = values[0][0]
            sem = 0
        else:
            order_param, _, sem = values[0]

        if 'M_G1_M' in atom_name:
            head_group['g1'].append((order_param, sem))
        elif 'M_G2_M' in atom_name:
            head_group['g2'].append((order_param, sem))
        elif 'M_G3_M' in atom_name:
            head_group['g3'].append((order_param, sem))
        elif 'M_G3N6' in atom_name:
            if 'C1' in atom_name:
                head_group['γ'].append((order_param, sem))
            elif 'C2' in atom_name:
                head_group['γ'].append((order_param, sem))
            elif 'C3' in atom_name:
                head_group['γ'].append((order_param, sem))
        elif 'M_G3C4' in atom_name:
            head_group['α'].append((order_param, sem))
        elif 'M_G3C5' in atom_name:
            head_group['β'].append((order_param, sem))

        if atom_name in sn1_keys:
            carbon_num = int(atom_name.split('C')[1].split('_')[0])
            sn1[carbon_num].append((order_param, sem))
        
        elif atom_name in sn2_keys:
            carbon_num = int(atom_name.split('C')[1].split('_')[0])
            sn2[carbon_num].append((order_param, sem))

    return head_group, sn1, sn2

def plot_data(data_files):
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    
    for json_file, color, label in data_files:
        json_data = load_json_data(json_file)
        head_group, sn1, sn2 = parse_data(json_data)

#headgroup data
        x_head_group = []
        y_head_group = []
        y_err_head_group = []

        for atom in ['γ', 'α', 'β', 'g3', 'g2', 'g1']:
            if atom in head_group and len(head_group[atom]) > 0:
                for order_param, sem in head_group[atom]:
                    x_head_group.append(atom)
                    y_head_group.append(order_param)
                    y_err_head_group.append(sem)

        if len(x_head_group) != len(y_head_group) or len(x_head_group) != len(y_err_head_group):
            print(f"Length mismatch in head group plot data: {len(x_head_group)}, {len(y_head_group)}, {len(y_err_head_group)}")
            continue

        ax1.errorbar(x_head_group, y_head_group, yerr=y_err_head_group, fmt='o', label=label, color=color) #, capsize=5
        ax1.set_title('Head Group')
        ax1.set_xlabel('Carbon')
        ax1.set_ylabel(r'$S_{CH}$', fontsize = 12)

#sn-1 data
        x_sn1 = []
        y_sn1 = []
        y_err_sn1 = []
        for carbon_num, hydrogens in sn1.items():
            for order_param, sem in hydrogens:
                x_sn1.append(carbon_num-1)
                y_sn1.append(abs(order_param))
                y_err_sn1.append(sem)
        ax2.errorbar(x_sn1, y_sn1, yerr=y_err_sn1, fmt='o-', label=label, color=color) #, capsize=5
        ax2.set_title('sn-1')
        ax2.set_xlabel('Carbon Number')
        ax2.set_ylabel(r'$|S_{CH}|$', fontsize = 12)

#sn-2 data
        x_sn2 = []
        y_sn2 = []
        y_err_sn2 = []
        for carbon_num, hydrogens in sn2.items():
            for order_param, sem in hydrogens:
                x_sn2.append(carbon_num-1)
                y_sn2.append(abs(order_param))
                y_err_sn2.append(sem)
        ax3.errorbar(x_sn2, y_sn2, yerr=y_err_sn2, fmt='o-', label=label, color=color) #, capsize=5
        ax3.set_title('sn-2')
        ax3.set_xlabel('Carbon Number')
        ax3.set_ylabel(r'$|S_{CH}|$', fontsize = 12)

    ax1.legend()
    ax2.legend()
    ax3.legend()
    fig.tight_layout()
    plt.savefig('ORDERPARAM_PLOT', dpi=300)
    plt.show()

def main():
    data_files = [
        ('Experiment_POPC_10.1039_c2cp42738a.json', 'black', 'Experiment'),
        # ('Experiment_POPC_nesterenko.json', 'black', 'Experiment'),

        # ('POPCOrderParameters_OpenFF_0.9_OG.json', 'blue', 'OpenFF'),
        # ('POPC_OpenFF_anneal.json', '#ADD8E6', 'OpenFF w/ Annealing'),

        # ('POPC_CHARMM_OG.json', 'red', 'CHARMM36'),
        # ('POPC_CHARMM_OG_ranonmysystem.json', 'orange', 'CHARMM36 Validation'),
        # ('POPC_CHARMM_control.json', '#FFE5B4', 'CHARMM36 Control'),                         ##FFBF00 AMBER

        ('POPC_Lipid17_OG.json', 'green', 'Lipid17'),
        # ('POPC_Lipid17_ranonmysystem.json', '#90EE90', 'Lipid17 Validation'),

        # ('POPC_Lipid21_anneal.json', 'purple', 'Lipid21 w/ Annealing'),
        # ('POPC_Lipid21_no_anneal.json', '#CBC3E3', 'Lipid21 w/o Annealing'),
        # ('POPC_Lipid21_BEST_unfinished.json', 'blue', 'Lipid21 Parmed'),
        ('POPC_Lipid21_newest.json', 'blue', 'Lipid21 in GROMACS'),

        
    ]
    plot_data(data_files)

if __name__ == "__main__":
    main()
