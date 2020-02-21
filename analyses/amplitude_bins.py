import numpy as np
import ROOT

from amplitudes.amplitude import Amplitude
from amplitudes.angular_dependencies import BelleS, BelleP, BelleD
from amplitudes.breit_wigner import BreitWigner

from analyses.collect_data import read_data_monte_carlo

from bins.utils import is_valid_dalitz_point
from bins.bins import create_bins


def plot_angular_function():
    m_dc = 1.86958
    m_pi = 0.13957
    m_kc = 0.493677

    mother_mass = m_dc
    fs_masses = [m_kc, m_pi, m_pi]

    bin_list, binning_x, binning_y = create_bins(m_dc, n_bins=200)

    # fcns = [osh, heb, hub, hob, uesen, halfBudalf]
    wave = BelleD(isobar_index=12, fs_masses=fs_masses)
    print("Amplitude created and functions assigned.")

    hist_real = ROOT.TH2D(f"Dalitz", f"Dalitz", len(binning_x) - 1, binning_x, len(binning_y) - 1, binning_y)

    for bin in bin_list:
        grid = bin.make_grid(mesh_width=0.01)
        grid = grid[is_valid_dalitz_point(grid, m_dc ** 2, m_kc ** 2, m_pi ** 2, m_pi ** 2)]
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
                        kin = np.array([np.power(mother_mass, 2), s_12, s_13])

                        val = wave.eval(kin)
                        # print(val)
                        hist_real.SetBinContent(i, j, val.real)

    hist_real.Draw('col')
    input()


def plot_breit_function():
    m_dc = 1.86958
    m_pi = 0.13957
    m_kc = 0.493677

    mother_mass = m_dc
    fs_masses = [m_kc, m_pi, m_pi]
    bachelor_mass = fs_masses[0]
    daughter_mass1 = fs_masses[1]
    daughter_mass2 = fs_masses[2]

    bin_list, (binning_x, binning_y) = create_bins(m_dc, n_bins=200)

    wave = BreitWigner(name='TestBreitWigner', mass=1, width=1, spin=0,
                       mother_mass=mother_mass,
                       bachelor_mass=bachelor_mass, daughter_mass1=daughter_mass1, daughter_mass2=daughter_mass2)

    print("Amplitude created and functions assigned.")

    hist_real = ROOT.TH2D(f"Dalitz", f"Dalitz", len(binning_x) - 1, binning_x, len(binning_y) - 1, binning_y)

    for bin in bin_list:
        grid = bin.make_grid(mesh_width=0.01)
        grid = grid[is_valid_dalitz_point(grid, m_dc ** 2, m_kc ** 2, m_pi ** 2, m_pi ** 2)]
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
                        val = wave.eval(s_13)
                        # print(val)
                        hist_real.SetBinContent(i, j, val.real)

    hist_real.Draw('col')
    input()


def plot_amplitude_bins():
    m_dc = 1.86958
    m_pi = 0.13957
    m_kc = 0.493677

    fs_masses = [m_kc, m_pi, m_pi]

    bin_list, binning_x, binning_y = create_bins(m_dc, n_bins=200)

    # fcns = [osh, heb, hub, hob, uesen, halfBudalf]
    amplitude = Amplitude(mother_mass=m_dc, fs_masses=fs_masses, mass=1, width=1)
    print("Amplitude created and functions assigned.")

    hist_real = ROOT.TH2D(f"Dalitz", f"Dalitz", len(binning_x) - 1, binning_x, len(binning_y) - 1, binning_y)
    hist_imag = ROOT.TH2D(f"Dalitz", f"Dalitz", len(binning_x) - 1, binning_x, len(binning_y) - 1, binning_y)

    for bin in bin_list:
        grid = bin.make_grid(mesh_width=0.001)
        grid = grid[is_valid_dalitz_point(grid, m_dc ** 2, m_kc ** 2, m_pi ** 2, m_pi ** 2)]
        # print(f'grid: {grid}')
        if len(grid) != 0:
            for grid_val in grid:
                # print(f'grid val: {grid_val}')
                borders = bin.get_borders()
                i_min = hist_real.GetXaxis().FindBin(borders[0])
                i_max = hist_real.GetXaxis().FindBin(borders[1])
                j_min = hist_real.GetYaxis().FindBin(borders[2])
                j_max = hist_real.GetYaxis().FindBin(borders[3])
                for i in range(i_min, i_max + 1):
                    for j in range(j_min, j_max + 1):
                        val = amplitude.eval(dalitz_point=grid_val)
                        # print(val)
                        hist_real.SetBinContent(i, j, val.real)
                        hist_imag.SetBinContent(i, j, val.imag)

    hist_real.Draw('col')
    input()

    # hist_imag.Draw('col')
    # input()


def plot_data_simple():
    m_dc = 1.86958
    m_pi = 0.13957
    m_kc = 0.493677

    data, data_size = read_data_monte_carlo()
    data = data.tolist()

    bin_list, binning_x, binning_y = create_bins(m_dc, n_bins=200)

    hist = ROOT.TH2D(f"Dalitz", f"Dalitz", len(binning_x) - 1, binning_x, len(binning_y) - 1, binning_y)

    progress_counter = 0
    for bin in bin_list:
        progress_counter += 1
        print(f'Progress: {int(100 * progress_counter / len(bin_list))}%')

        for point in data:
            if bin.contains(*point):
                bin.increment_value()

        borders = bin.get_borders()
        i_min = hist.GetXaxis().FindBin(borders[0])
        i_max = hist.GetXaxis().FindBin(borders[1])
        j_min = hist.GetYaxis().FindBin(borders[2])
        j_max = hist.GetYaxis().FindBin(borders[3])
        for i in range(i_min, i_max + 1):
            for j in range(j_min, j_max + 1):
                hist.SetBinContent(i, j, bin.value)

    hist.Draw('col')
    input()


if __name__ == "__main__":
    plot_breit_function()
    # plot_angular_function()
    # plot_data_simple()

    # plot_amplitude_bins()

