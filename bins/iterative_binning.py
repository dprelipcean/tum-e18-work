import numpy as np
import pickle

from amplitudes.constants import m_dc, m_kc, m_pi, fs_masses
from analyses.collect_data import read_data_monte_carlo
from analyses.amplitude_bins import plot_values
from bins.bins_main import RectangularBin
from bins.utils import is_valid_dalitz_point_pointwise


def create_bins_matrix(m_dc, n_bins):
    """Create the bins for tor the Dalitz plot."""
    binning_x = np.linspace(0., m_dc ** 2, n_bins)
    binning_y = np.linspace(0., m_dc ** 2, n_bins)

    bin_matrix = []
    for i in range(len(binning_x)-1):
        bin_column = []
        for j in range(len(binning_y)-1):
            mybin = RectangularBin(binning_x[i], binning_x[i + 1], binning_y[j], binning_y[j + 1])

            if is_valid_dalitz_point_pointwise(mybin.position, m_dc ** 2, fs_masses[0] ** 2, fs_masses[1] ** 2, fs_masses[2] ** 2):
                bin_column.append(mybin)

        bin_matrix.append(bin_column)

    return bin_matrix

n_bins = 300
bin_matrix = create_bins_matrix(m_dc, n_bins=n_bins)


def get_bin_value(s_12, s_13, bin):
    return bin.value


def iteration_method(my_bin_matrix, bin_threshold_value):
    i = len(my_bin_matrix) - 3
    while i >= 0:
        j = min(len(my_bin_matrix[i]) - 3, len(my_bin_matrix[i + 1]) - 2)
        while j >= 0:

            bin = my_bin_matrix[i][j]
            bin_upper = my_bin_matrix[i][j + 1]
            bin_right = my_bin_matrix[i + 1][j]
            bin_diago = my_bin_matrix[i + 1][j + 1]

            bin_sum_values = bin.value + bin_upper.value + bin_right.value + bin_diago.value
            if bin_sum_values < bin_threshold_value:
                bin.redefine_borders(bin.get_borders()[0], bin_diago.get_borders()[1],
                                     bin.get_borders()[2], bin_diago.get_borders()[3])
                bin.value = bin_sum_values

                del my_bin_matrix[i + 1][j + 1]
                del my_bin_matrix[i + 1][j]
                del my_bin_matrix[i][j + 1]

            j -= 2
        i -= 2


def iteration_method_2(my_bin_matrix, bin_threshold_value):
    i = len(my_bin_matrix) - 2
    while i >= 0:
        j = min(len(my_bin_matrix[i]) - 2, len(my_bin_matrix[i + 1]) - 2)
        while j >= 0:

            bin = my_bin_matrix[i][j]
            bin_upper = my_bin_matrix[i][j + 1]

            bin_sum_values = bin.value + bin_upper.value

            if bin_sum_values < bin_threshold_value:
                bin.redefine_borders(bin.get_borders()[0], bin_upper.get_borders()[1],
                                     bin.get_borders()[2], bin_upper.get_borders()[3])
                bin.value = bin_sum_values

                del my_bin_matrix[i][j + 1]

            j -= 1
        i -= 1

    j = min(len(my_bin_matrix[i]) - 2, len(my_bin_matrix[i + 1]) - 2)
    while j >= 0:
        i = len(my_bin_matrix) - 2
        while i >= 0:

            bin = my_bin_matrix[i][j]
            bin_right = my_bin_matrix[i + 1][j]

            bin_sum_values = bin.value + bin_right.value
            if bin_sum_values < bin_threshold_value:
                bin.redefine_borders(bin.get_borders()[0], bin_right.get_borders()[1],
                                     bin.get_borders()[2], bin_right.get_borders()[3])
                bin.value = bin_sum_values

                del my_bin_matrix[i + 1][j]
            j -= 1
        i -= 1

def plot_data_iteratively():
    bin_threshold_value = 10
    max_threshold_value = 200
    threshold_increment = 10
    flatten = lambda l: [item for sublist in l for item in sublist]


    print(f"Reading data.")

    data, data_size = read_data_monte_carlo()
    data = data.tolist()

    print(f"Computing data values in plot.")
    bins_list_from_matrix = flatten(bin_matrix)
    print(len(bins_list_from_matrix))

    for bin_column in bin_matrix:
        for bin in bin_column:
            for point in data:
                if bin.contains(*point):
                    bin.increment_value()

    for iteration_step in range(int(max_threshold_value / threshold_increment)):
        print(f"Computing iteration step: {iteration_step}")
        if bin_threshold_value < max_threshold_value:
            bin_threshold_value += threshold_increment

        iteration_method_2(bin_matrix, bin_threshold_value)

        bins_list_from_matrix = flatten(bin_matrix)
        print(len(bins_list_from_matrix))

    with open('./binning_list.bin', 'wb') as binning_file:
        pickle.dump(bins_list_from_matrix, binning_file)

    plot_values(value_function=get_bin_value, bins=bins_list_from_matrix)


if __name__ == '__main__':
    plot_data_iteratively()
