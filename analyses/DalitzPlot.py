# -*- coding: utf-8 -*-
#
# Copyright (C) 2019, 2020 TUM E18 Research Group.

"""Dalitz plotting."""

import numpy as np

from amplitudes.breit_wigner import BreitWigner

from analyses.collect_data import read_data_monte_carlo
from bins.bins_main import create_bins
from bins.chi2model import DalitzChi2model
from bins.utils import is_valid_dalitz_point, hub
from amplitudes.amplitude import Amplitude
from amplitudes.angular_dependencies import BelleS, BelleP, BelleD

from data.metadata import DATA_FILE_CHANNELS

from amplitudes.constants import m_dc, fs_masses


def initialize_functions():
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

    return fcns


def main():

    bin_list = create_bins(m_dc, n_bins=200)[0]

    fcns = initialize_functions()

    print("Amplitude created and functions assigned.")

    c2model = DalitzChi2model(m_dc, fs_masses, bin_list, 0.01, fcns, 10)
    print("Chi2 model created.")

    c2model.make_valid_bins()
    c2model.make_integrals()

    # m2s = generate_random_data(m_pi, m_kc, m_dc,< n_data=10000)
    m2s, _ = read_data_monte_carlo()

    valid = is_valid_dalitz_point(m2s, m_dc ** 2, fs_masses[0] ** 2, fs_masses[1] ** 2, fs_masses[2] ** 2)
    print(f"nGOOD ={np.sum(valid)}")

    m2s = m2s[valid, :]
    weights = 1. + hub(m2s) ** 2

    mxx = np.max(weights)

    weights -= np.random.random(weights.shape) * mxx

    m2s = m2s[weights > 0., :]

    print(f"nDeWeight ={np.sum(weights > 0)}")

    c2model.load_data(m2s[:, 0], m2s[:, 1])

    # params = np.random.uniform(-1.,1.,2*len(fcns))
    params = np.zeros(2 * len(fcns))

    from iminuit import Minuit

    m = Minuit.from_array_func(c2model.eval, params, error=0.5)
    m.migrad()

    vals = np.array([m.values['x' + str(i)] for i in range(2 * len(fcns))])

    h1 = c2model.make_theo_hist(vals)
    h1.Draw('col')
    input()
    h2 = c2model.make_data_hist()
    h2.Draw('col')
    input()
    print(vals)

    ntfrr = ((vals[0] + 1.j * vals[1]) * (vals[2] - 1.j * vals[3])) ** 2
    print(f"{ntfrr / abs(ntfrr)} should be real (phase of pi/2)")


if __name__ == "__main__":
    main()
