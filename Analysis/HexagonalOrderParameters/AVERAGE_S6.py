
# WIP

import numpy as np

def read_data_from_file(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip():
                values = line.split()
                index = int(values[0])
                second_column = float(values[1])
                third_column = float(values[2])
                data.append((index, second_column, third_column))
    return data

def calculate_averages_and_sem(file_path):
    data = read_data_from_file(file_path)

    second_column = [row[1] for row in data]
    third_column = [row[2] for row in data]

    average_second_column = np.mean(second_column)
    average_third_column = np.mean(third_column)

    std_second_column = np.std(second_column, ddof=1)
    std_third_column = np.std(third_column, ddof=1)

    sem_second_column = std_second_column / np.sqrt(len(second_column))
    sem_third_column = std_third_column / np.sqrt(len(third_column))

    return (average_second_column, sem_second_column), (average_third_column, sem_third_column)

file_path = "s6.xvg"

(average_second_column, sem_second_column), (average_third_column, sem_third_column) = calculate_averages_and_sem(file_path)

print(f"Average S6: {average_second_column:.4f} +/- {sem_second_column:.4f}")
print(f"Average computed fraction of chains with S6 > 0.72: {average_third_column:.4f} +/- {sem_third_column:.4f}")