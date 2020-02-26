# -*- coding: utf-8 -*-
#
# Copyright (C) 2019, 2020 TUM E18 Research Group.

"""Classes implementing the Breit Wigner function.

The Breit-Wigner is always called, with the invariant mass-square of the
two-particle combination, that forms the Isobar. The other combinations
have to be chosen accordingly in the angular part.
"""
import numpy as np

from amplitudes.plotting_scripts import plot


class BreitWigner:
    """Implementation of the BreitWigner function."""

    def __init__(self, name, mass, width, spin, mother_mass, bachelor_mass, daughter_mass1, daughter_mass2,
                 rr=1.5, rd=5.0):
        """

        Parameters
        ----------
        name
        mass
        width
        spin
        mother_mass
        bachelor_mass
        daughter_mass1
        daughter_mass2
        rr: float, optional
            Value taken from the C++ equivalent code.
            Defaults to 1.5.
        rd: float, optional
            Value taken from the C++ equivalent code.
            Defaults to 5.
        """
        self.name = name

        self._check_conditions(mother_mass, daughter_mass1, daughter_mass2, spin)

        self.mass = mass
        self.width = width
        self.spin = spin

        self.daughter_mass1 = daughter_mass1
        self.daughter_mass2 = daughter_mass2

        self.bachelor_mass = bachelor_mass
        self.mother_mass = mother_mass

        self.rr = rr
        self.rd = rd

    def _check_conditions(self, mass, daughter_mass1, daughter_mass2, spin):
        """Check conditions on input data."""
        self._check_condition_daughter_masses(mass, daughter_mass1, daughter_mass2)
        self._check_condition_spin(spin)

    def _check_condition_daughter_masses(self, mass, daughter_mass1, daughter_mass2):
        """Check the physical condition that the mass of the two daughter cannot be greater than the initial mass."""
        # print(f'Daughter mass1: {daughter_mass1}\n'
        #       f'Daughter mass2: {daughter_mass2}\n'
        #       f'Mass: {mass}\n')
        if (daughter_mass1 + daughter_mass2) > mass:
            raise Exception(f"BELLEbreitWigner::BELLEbreitWigner(...): ERROR: On shell resonance mass of {self.name} "
                            f"too light for decay into daughter particles:{mass} < {daughter_mass1}+"
                            f"{daughter_mass2}")

    @staticmethod
    def _check_condition_spin(spin):
        """Check the spin for implemented scenarios."""
        if spin > 2:
            raise Exception(f"BELLEbreitWigner::BELLEbreitWigner(...): ERROR: Spin > 2 not supported yet")

    def eval(self, s12, constant_gamma=True):
        """Evaluate the Breit-Wigner function."""
        m12 = np.power(s12, 1/2)
        # self._check_condition_daughter_masses(m12, self.daughter_mass1, self.daughter_mass2)

        # ToDo: check why is this needed
        # if (m12 > self.mother_mass - self.bachelor_mass) or (m12 < self.daughter_mass2 + self.daughter_mass2):
        #     print("Masses condition not satisfied.")
        #     return 0 + 0 * 1j

        s = np.power(self.mother_mass, 2)
        sr = np.power(self.mass, 2)
        s_daughter1 = np.power(self.daughter_mass1, 2)
        s_daughter2 = np.power(self.daughter_mass2, 2)
        s_batch = np.power(self.bachelor_mass, 2)

        val_1 = np.power(sr - s_daughter1 - s_daughter2, 2) - 4 * s_daughter1 * s_daughter2
        val_2 = np.power(s12 - s_daughter1 - s_daughter2, 2) - 4 * s_daughter1 * s_daughter2

        p_r = np.power(val_1, 0.5) / 2 / self.mass
        p_ab = np.power(val_2, .5) / 2 / m12

        # p_d = np.power(np.power(s - sr - s_batch, 2) - 4 * sr * s_batch, .5) / 2 / self.mother_mass
        s_val = np.power(s - s12 - s_batch, 2) - 4 * s12 * s_batch
        # Do (sign * abs) to avoid RuntimeWarning: invalid value encountered in power
        p_abc = np.sign(s_val) * np.power(np.abs(s_val), .5) / 2 / self.mother_mass

        if self.spin == 0:
            fr = 1.
            fd = 1.
        elif self.spin == 1:
            fr = np.power((np.power(self.rr * p_r, 2) + 1) / (np.power(self.rr * p_ab, 2) + 1), .5)
            # fd = np.power((np.power(self.RD * p_d, 2)+1) / (np.power(self.RD * p_abc, 2)+1), .5)
            fd = np.power(1. / (np.power(self.rd * p_abc, 2) + 1), .5)
        else:  # self.spin == 2:
            xr = self.rr * self.rr * p_r * p_r
            x_ab = self.rr * self.rr * p_ab * p_ab

            fr = np.power((np.power(xr-3., 2) + 9 * xr) / (np.power(x_ab-3., 2) + 9 * x_ab), .5)

            # xD   = self.RD * self.RD * p_d * p_d
            x_abc = np.power(self.rd, 2) * np.power(p_abc, 2)
            # fd = np.power((np.power(xD-3., 2) + 9 * xD) / (np.power(x_abc-3., 2) + 9 * x_abc), .5)
            fd = np.power(1. / (np.power(x_abc-3., 2) + 9 * x_abc), .5)

        gamma = self.width * self.mass / m12 * np.power(fr, 2) * np.power(p_ab / p_r, 2 * self.spin + 1)

        # Use this for complex BreitWigner:
        ret_val = complex(fr * fd, 0.) / complex(np.power(self.mass, 2) - s12, - self.mass * gamma)

        # Use this for real BreitWigner:
        # ret_val = np.abs(complex(fr * fd, 0.) / ((np.power(self.mass, 2) - s12) + np.power(self.mass * gamma, 2)))

        # if not ret_val.real or not ret_val.imag:
        #     print(f"s12 = {s12} sDaughter1 = {s_daughter1} sDaughter2 = {s_daughter2} p_r = {p_r} p_ab = {p_ab} "
        #           f"fr = {fr} fd = {fd} gamma = {gamma}")
        # print(f'ret_val: {ret_val}')

        return ret_val


def plot_breit_wigner():
    """Plot the BreitWigner function."""


    m_dc = 1.86958
    m_pi = 0.13957
    m_kc = 0.493677

    fs_masses = [m_kc, m_pi, m_pi]

    bachelor_mass = fs_masses[0]
    daughter_mass1 = fs_masses[1]
    daughter_mass2 = fs_masses[2]


    test_points = np.linspace(0, 4, 200)
    val_points = list()

    breit_wigner_points_real = list()
    breit_wigner_points_imaginary = list()

    breit_wigner_points_amplitude = list()

    test_mass = 1

    for val in test_points:
        breit_wigner = BreitWigner(name='TestBreitWigner', mass=test_mass, width=0.1, spin=0,
                                   mother_mass=m_dc,
                                   bachelor_mass=bachelor_mass,
                                   daughter_mass1=daughter_mass1, daughter_mass2=daughter_mass2)
        ret_val = breit_wigner.eval(val)

        breit_wigner_points_real.append(ret_val.real)
        breit_wigner_points_imaginary.append(ret_val.imag)

        breit_wigner_points_amplitude.append(abs(ret_val))

        val_points.append(val)

    plot(val_points, breit_wigner_points_real, breit_wigner_points_imaginary, breit_wigner_points_amplitude)


if __name__ == '__main__':
    plot_breit_wigner()
