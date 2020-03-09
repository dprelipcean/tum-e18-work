import numpy as np
import ROOT

from amplitudes.amplitude import Amplitude
from amplitudes.angular_dependencies import BelleS, BelleP, BelleD
from amplitudes.breit_wigner import BreitWigner

from analyses.collect_data import read_data_monte_carlo
from amplitudes.constants import m_dc, m_kc, m_pi, fs_masses

from bins.utils import is_valid_dalitz_point
from bins.bins_main import create_bins
from bins.iterative_binning import create_bins_matrix

mother_mass = m_dc

bachelor_mass = fs_masses[0]
daughter_mass1 = fs_masses[1]
daughter_mass2 = fs_masses[2]

n_bins = 100
bins_list, (binning_x, binning_y) = create_bins(m_dc, n_bins=n_bins)
bin_matrix = create_bins_matrix(m_dc, n_bins=n_bins)


def plot_values(wave=None, value_function=None, bins=None):

    if bins:
        bin_list = bins
    else:
        bin_list = bins_list

    hist_real = ROOT.TH2D(f"Dalitz", f"Dalitz", len(binning_x) - 1, binning_x, len(binning_y) - 1, binning_y)

    for bin in bin_list:
        grid = bin.make_grid(mesh_width=0.01)
        grid = grid[is_valid_dalitz_point(grid, m_dc ** 2, m_kc ** 2, m_pi ** 2, m_pi ** 2)]

        if bins:
            wave = bin

        # print(f'grid: {grid}')
        if len(grid) != 0:
            for grid_val in grid:
                # print(f'grid val: {grid_val}')
                s_12 = grid_val[0]
                s_13 = grid_val[1]

                borders = bin.get_borders()
                i_min = hist_real.GetXaxis().FindBin(borders[0])
                i_max = hist_real.GetXaxis().FindBin(borders[1])
                j_min = hist_real.GetYaxis().FindBin(borders[2])
                j_max = hist_real.GetYaxis().FindBin(borders[3])
                for i in range(i_min, i_max + 1):
                    for j in range(j_min, j_max + 1):
                        val = value_function(s_12, s_13, wave)
                        hist_real.SetBinContent(i, j, val.real)

    hist_real.Draw('col')
    input()


def get_angular_function_value(s_12, s_13, wave):
    kin = np.array([np.power(mother_mass, 2), s_12, s_13])
    val = wave.eval(kin)
    return val


def plot_angular_function():
    wave = BelleD(isobar_index=12, fs_masses=fs_masses)
    plot_values(wave, get_angular_function_value)


def get_breit_value(s_12, s_13, wave):
    val = wave.eval(s_13)
    return val


def plot_breit_wigner():
    wave = BreitWigner(name='TestBreitWigner', mass=1, width=1, spin=0,
                       mother_mass=mother_mass,
                       bachelor_mass=bachelor_mass, daughter_mass1=daughter_mass1, daughter_mass2=daughter_mass2)
    plot_values(wave, get_breit_value)


def get_amplitude_value(s_12, s_13, amplitude):
    grid_val = [s_12, s_13]
    val = amplitude.eval(dalitz_point=grid_val)
    return val.real


def plot_amplitude():
    amplitude = Amplitude(resonance_model=BreitWigner, angular_dependence_class=BelleS,
                          mother_mass=m_dc, fs_masses=fs_masses, mass=1, width=1)
    plot_values(amplitude, get_amplitude_value)


def get_bin_value(s_12, s_13, bin):
    return bin.value


def plot_data_simple():
    data, data_size = read_data_monte_carlo(data_range=10000)
    data = data.tolist()

    for bin in bins_list:
        for point in data:
            if bin.contains(*point):
                bin.increment_value()

    plot_values(value_function=get_bin_value, bins=bins_list)


def plot_data_iteratively():
    print(f"Reading data.")

    data, data_size = read_data_monte_carlo(data_range=None)
    data = data.tolist()

    print(f"Computing data values in plot.")
    for bin_column in bin_matrix:
        for bin in bin_column:
            for point in data:
                if bin.contains(*point):
                    bin.increment_value()

    # flatten = lambda l: [item for sublist in l for item in sublist]
    # plot_values(value_function=get_bin_value, bins=flatten(bin_matrix))

    for iteration_step in range(20):
        print(f"Computing iteration step: {iteration_step}")
        i = len(bin_matrix) - 3
        bin_threshold_value = 10
        while i >= 0:
            j = min(len(bin_matrix[i]) - 3, len(bin_matrix[i+1]) -2)
            while j >= 0:
                # print(f'{i}, {j}')
                # print(f'lengths: {len(bin_matrix)} {len(bin_matrix[i])} {len(bin_matrix[i+1])}')

                bin = bin_matrix[i][j]
                bin_upper = bin_matrix[i][j + 1]
                bin_right = bin_matrix[i + 1][j]
                bin_diago = bin_matrix[i + 1][j + 1]

                bin_sum_values = bin.value + bin_upper.value + bin_right.value + bin_diago.value
                if bin_sum_values < bin_threshold_value:
                    bin.borders = [bin.get_borders()[0], bin_diago.get_borders()[1], bin.get_borders()[2],
                                   bin_diago.get_borders()[3]]
                    bin.value = bin_sum_values

                    del bin_matrix[i + 1][j + 1]
                    del bin_matrix[i + 1][j]
                    del bin_matrix[i][j + 1]

                j -= 2
            i -= 2
        

    flatten = lambda l: [item for sublist in l for item in sublist]
    plot_values(value_function=get_bin_value, bins=flatten(bin_matrix))


if __name__ == "__main__":
    # plot_breit_function()
    # plot_angular_function()
    # plot_data_simple()
    # plot_amplitude_bins()
    plot_data_iteratively()

