# -*- coding: utf-8 -*-
#
# Copyright (C) 2019, 2020 TUM E18 Research Group.

"""Dalitz plotting."""

import numpy as np

from amplitudes.amplitude import Amplitude
from amplitudes.angular_dependencies import BelleS, BelleP, BelleD
from amplitudes.breit_wigner import BreitWigner
from amplitudes.constants import m_dc, fs_masses

from analyses.collect_data import read_data_monte_carlo
from analyses.perform_fit import perform_fit

from bins.bins_main import define_binning
from bins.chi2model import DalitzChi2model
from bins.utils import is_valid_dalitz_point, compute_weights

from data.metadata import DATA_FILE_CHANNELS



def initialize_functions():
    """Initialize the amplitude functions for the resonances."""
    fcns = list()
    for wave in DATA_FILE_CHANNELS.keys():
        if wave == 's_wave':
            angular_dependence_class = BelleS
        elif wave == 'p_wave':
            angular_dependence_class = BelleP
        else:
            angular_dependence_class = BelleD

        for resonance in DATA_FILE_CHANNELS[wave].values():
            mass = resonance['mass']
            width = resonance['width']

            amplitude = Amplitude(mass=mass, width=width, mother_mass=m_dc, fs_masses=fs_masses,
                                  resonance_model=BreitWigner,
                                  angular_dependence_class=angular_dependence_class)
            fcns.append(amplitude.eval)
    print("Amplitude created and functions assigned.")

    return fcns


def main(plot_dalitz_plane=True):

    bin_list = define_binning()
    fcns = initialize_functions()

    chi2model = DalitzChi2model(m_dc, fs_masses, bin_list, 0.01, fcns, 10)
    print("Chi2 model created.")

    chi2model.make_valid_bins()
    chi2model.make_integrals()

    m2s, _ = read_data_monte_carlo()

    valid = is_valid_dalitz_point(m2s, m_dc ** 2, fs_masses[0] ** 2, fs_masses[1] ** 2, fs_masses[2] ** 2)
    print(f"nGOOD ={np.sum(valid)}")

    m2s = m2s[valid, :]

    weights = compute_weights(m2s)

    m2s = m2s[weights > 0., :]

    print(f"nDeWeight ={np.sum(weights > 0)}")

    chi2model.load_data(m2s[:, 0], m2s[:, 1])
    print(f"Data loaded")

    vals = perform_fit(chi2model, fcns)

    if plot_dalitz_plane:
        h1 = chi2model.make_theo_hist(vals)
        h1.Draw('col')
        input()
        # h2 = chi2model.make_data_hist()
        # h2.Draw('col')
        # input()


if __name__ == "__main__":
    main()
