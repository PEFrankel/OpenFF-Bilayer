import json
import matplotlib.pyplot as plt
import numpy as np

def read_json_to_dict(filename):
    with open(filename, 'r') as file:
        data_dict = json.load(file)
    return data_dict

def filter_data_points(data_dict, interval=200): # plot data for every 200 ps since runs saved data differently
    times = list(map(float, data_dict.keys()))
    apl_values = list(data_dict.values())
    
    filtered_times = [times[i] for i in range(0, len(times), interval // 10)]
    filtered_apl_values = [apl_values[i] for i in range(0, len(apl_values), interval // 10)]
    
    return filtered_times, filtered_apl_values

def plot_graph(file_list, color_map, legend_map):

    # plt.figure(figsize=(10, 6))
    
    for file in file_list:
        data_dict = read_json_to_dict(file)
        times, apl_values = filter_data_points(data_dict)
        
        color = color_map.get(file, 'black')
        legend_label = legend_map.get(file, file)
        plt.plot(times, apl_values, color=color, alpha=0.5, label=legend_label)

        coeffs = np.polyfit(times, apl_values, 1)
        trendline = np.polyval(coeffs, times)
        plt.plot(times, trendline, color=color, linestyle='--')
    
    plt.xlabel("Time (ps)")
    plt.ylabel("APL (Å²)")
    plt.xlim([0,500000])
    plt.title('APL vs Time')
    plt.legend()

    plt.savefig('APL_all2.png', dpi = 300, bbox_inches='tight')
    plt.show()

#  main()
file_list = ['apl_packmol.json','apl_macrog.json','apl_slipids.json']
color_map = {
    'apl_packmol.json': '#0333b0',
    'apl_macrog.json': 'green',
    'apl_slipids.json': '#ee7f17',
}
legend_map = {
    'apl_packmol.json': 'OpenFF (Eq time = 4.05)',
    'apl_macrog.json': 'MacRog (Eq time = 0.32)',
    'apl_slipids.json': 'Slipids (Eq time = 0.09)'
}
plot_graph(file_list, color_map, legend_map)
