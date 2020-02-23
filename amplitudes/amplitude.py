# -*- coding: utf-8 -*-
#
# Copyright (C) 2019, 2020 TUM E18 Research Group.

"""Classes implementing the resonance amplitudes.

The variable naming always refers to a decay of a mother M into three final-state
particles 1,2,3. So mMother is the mass of the decaying particle and
fs_masses is a vector containing the three masses of the final state
particles {m1,m2,m3}.

Now we model the process not as direct decay [M -> 1 + 2 + 3], but
through two two-particle decays: [M -> Bachelor + Isobar] and [Isobar ->
daughter1 + daughter2], where the Isobar (an intermediate 2-body
resonance) can either be formed by the final-state particles (12), (13),
or (23). Therefore the daughter masses and the bachelor mass are all
elements of fs_masses, but the assignment depends on which two
particles form the Isobar.
"""

import numpy as np


class Amplitude:
    """Class implementing the amplitude as a combination of Breit Wigner and angular functions."""

    def __init__(self, mother_mass, fs_masses, mass, width, resonance_model=None, angular_dependence_class=None,
                 identical_particles=True):

        self.identical_particles = identical_particles

        self.mother_mass = mother_mass

        self.fs_masses = fs_masses

        self.mass = mass
        self.width = width

        self.angular_dependence_class = angular_dependence_class
        self.resonance_model = resonance_model

        name = 'Test'
        spin = 0
        rr = 0
        rd = 0

        bachelor_mass = self.fs_masses[0]
        daughter_mass1 = self.fs_masses[1]
        daughter_mass2 = self.fs_masses[2]

        self.breit_wigner_12 = self.resonance_model(name, self.mass, self.width, spin,
                                                    self.mother_mass, bachelor_mass,
                                                    daughter_mass1, daughter_mass2, rr, rd)

        self.partial_wave_12 = self.angular_dependence_class(isobar_index=12, fs_masses=self.fs_masses)

        self.breit_wigner_13 = self.resonance_model(name, self.mass, self.width, spin, self.mother_mass, bachelor_mass,
                                                    daughter_mass2, daughter_mass1, rr, rd)

        self.partial_wave_13 = self.angular_dependence_class(isobar_index=13, fs_masses=self.fs_masses)

    def eval(self, grid=None, dalitz_point=None):
        """Evaluate the function.

        Parameters
        ----------
        grid: ndarray
            Array containing Dalitz points as pairs, e.g. (m12, m23)
        dalitz_point: tuple

        Returns
        -------
        out: ndarray
            Array containing the evaluated output from the input Dalitz points grid.
        """

        if grid is not None:
            return_value = list()

            for dalitz_point in grid:
                val = self.eval(dalitz_point=dalitz_point)
                return_value.append(val)

            output = np.asarray(return_value)

        elif dalitz_point is not None:
            s_12 = dalitz_point[0]
            s_13 = dalitz_point[1]

            kin = np.array([np.power(self.mother_mass, 2), s_12, s_13])

            val = self._compute_value(s_12, s_13, kin)

            output = val

        else:
            output = None

        return output

    def _compute_value(self, s_12, s_13, kin):
        """Do the explicit computation."""
        if self.identical_particles:
            val = self.breit_wigner_12.eval(s_12) * self.partial_wave_12.eval(kin) \
                  + self.breit_wigner_13.eval(s_13) * self.partial_wave_13.eval(kin)
        else:
            val = None

        return val
