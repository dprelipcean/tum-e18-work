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

        self._check_conditions(mass, daughter_mass1, daughter_mass2, spin)

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
        self._check_condition_daughter_masses(daughter_mass1, daughter_mass2, mass)
        self._check_condition_spin(spin)

    def _check_condition_daughter_masses(self, daughter_mass1, daughter_mass2, mass):
        """Check the physical condition that the mass of the two daughter cannot be greater than the initial mass."""
        # ToDo
        # if (daughter_mass1 + daughter_mass2) > mass:
        #     raise Exception(f"BELLEbreitWigner::BELLEbreitWigner(...): ERROR: On shell resonance mass of {self.name} "
        #                     f"too light for decay into daughter particles:{mass} < {daughter_mass1}+"
        #                     f"{daughter_mass2}")
        pass

    @staticmethod
    def _check_condition_spin(spin):
        """Check the spin for implemented scenarios."""
        if spin > 2:
            raise Exception(f"BELLEbreitWigner::BELLEbreitWigner(...): ERROR: Spin > 2 not supported yet")

    def eval(self, s12):
        """Evaluate the Breit-Wigner function."""
        m12 = np.power(s12, 2)

        self._check_condition_daughter_masses(self.daughter_mass1, self.daughter_mass2, m12)

        if m12 > self.mother_mass - self.bachelor_mass:
            # ToDo: raise error or return mock value
            # raise Exception
            return 1 + 1 * 1j

        s = np.power(self.mother_mass, 2)
        sr = np.power(self.mass, 2)
        s_daughter1 = np.power(self.daughter_mass1, 2)
        s_daughter2 = np.power(self.daughter_mass2, 2)
        s_batch = np.power(self.bachelor_mass, 2)

        p_r = np.power(np.power(sr - s_daughter1 - s_daughter2, 2)
                       - 4 * s_daughter1 * s_daughter2, .5) / 2 / self.mass
        p_ab = np.power(np.power(s12 - s_daughter1 - s_daughter2, 2) - 4 * s_daughter1 * s_daughter2, .5) / 2 / m12

        # p_d = np.power(np.power(s - sr - s_batch, 2) - 4 * sr * s_batch, .5) / 2 / self.mother_mass
        p_abc = np.power(np.power(s - s12 - s_batch, 2) - 4 * s12 * s_batch, .5) / 2 / self.mother_mass

        fr = 1.
        fd = 1.

        if self.spin == 1:
            fr = np.power((np.power(self.rr * p_r, 2) + 1) / (np.power(self.rr * p_ab, 2) + 1), .5)
            # fd = np.power((np.power(self.RD * p_d, 2)+1) / (np.power(self.RD * p_abc, 2)+1), .5)
            fd = np.power(1. / (np.power(self.rd * p_abc, 2) + 1), .5)
        elif self.spin == 2:
            xr = self.rr * self.rr * p_r * p_r
            x_ab = self.rr * self.rr * p_ab * p_ab

            fr = np.power((np.power(xr-3., 2) + 9 * xr) / (np.power(x_ab-3., 2) + 9 * x_ab), .5)

            # xD   = self.RD * self.RD * p_d * p_d
            x_abc = np.power(self.rd, 2) * np.power(p_abc, 2)
            # fd = np.power((np.power(xD-3., 2) + 9 * xD) / (np.power(x_abc-3., 2) + 9 * x_abc), .5)
            fd = np.power(1. / (np.power(x_abc-3., 2) + 9 * x_abc), .5)

        gamma = self.width * self.mass / m12 * fr * fr * np.power(p_ab / p_r, 2 * self.spin + 1)

        ret_val = complex(fr * fd, 0.) / complex(np.power(self.mass, 2) - s12, - self.mass * gamma)

        if not ret_val.real or not ret_val.imag:
            print(f"s12 = {s12} sDaughter1 = {s_daughter1} sDaughter2 = {s_daughter2} p_r = {p_r} p_ab = {p_ab} "
                  f"fr = {fr} fd = {fd} gamma = {gamma}")

        return ret_val


class KinematicSignature:
    """Class for checking the kinematic signature."""
    nKin = 3


class AngularDependence:
    """Parent class for angular dependence functions."""

    valid_isobar_indexes = [12, 13, 23]

    def __init__(self):
        pass

    def _check_condition_isobarindex(self, isobar_index):
        if isobar_index not in self.valid_isobar_indexes:
            raise Exception(f"None of the three possible isobar masses match (12, 13 and 23) the given value: "
                            f"{isobar_index}")

    @staticmethod
    def _check_condition_fsmasses(fs_masses):
        if len(fs_masses) != 3:
            raise Exception(f"Number of final state masses given is not three.")

    @staticmethod
    def _check_condition_kinematic(kin):
        """Check the kinematic signature condition."""
        if len(kin) != KinematicSignature.nKin:
            raise Exception(f"BelleS::eval(...) Number of kinematic variables does not match ({kin.size}) + "
                            f"{KinematicSignature.nKin}")
        else:
            return True


class BelleAngDep(AngularDependence):
    """Class for implementing the Belle AngularDependence functions."""

    def __init__(self, isobar_index, fs_masses):
        """

        Parameters
        ----------
        isobar_index: 12, 13, 23
            Integer index denoting how is the isobar formed from the decaying particles.
        fs_masses: ndarray
            Vector containing the three masses of the final state particles (m1, m2, m3)
        """
        super(BelleAngDep, self).__init__()
        self._check_condition_isobarindex(isobar_index)
        self._check_condition_fsmasses(fs_masses)

        self.isobar_index = isobar_index
        self.fs_masses = fs_masses

    def set_fs_masses(self, new_masses):
        """Set new masses for the decaying particles.

        Parameters
        ----------
        new_masses: ndarray
        """
        if new_masses.size != 3:
            raise Exception(f"BelleP::setFSmasses(...)", "Number of masses given is not three: {newMasses.size}")
        self.fs_masses = new_masses

    def compute_s_kpi_wrong(self, m_d2, kin):
        """Compute the

        Parameters
        ----------
        m_d2: float

        kin

        Returns
        -------

        """
        return m_d2 + np.power(self.fs_masses[0], 2) + np.power(self.fs_masses[1], 2) \
                    + np.power(self.fs_masses[2], 2) - kin[1] - kin[2]

    def compute_masses(self, kin, s_k_pi_wrong):
        """

        Parameters
        ----------
        kin: ndarray
            The convention for the argument "kin" is: {mMother^2, m_12^2, m_13^2}
            so for an isobar formed by (12), mAB2 is set to kin[1] and so on.
        s_k_pi_wrong:
            # Todo: ?

        Returns
        -------
        m_a2: float
        m_b2:float
        m_c2: float
        m_ab2: float
        m_ac2: float
        m_bc2:float
        """
        if self.isobar_index == 12:  # ((piRight, Ks), piWrong)
            m_a2 = np.power(self.fs_masses[0], 2)
            m_b2 = np.power(self.fs_masses[1], 2)
            m_c2 = np.power(self.fs_masses[2], 2)

            m_ab2 = kin[1]
            m_ac2 = kin[2]
            m_bc2 = s_k_pi_wrong
        elif self.isobar_index == 13:  # ((piRight, piWrong), Ks)
            m_a2 = np.power(self.fs_masses[2], 2)
            m_b2 = np.power(self.fs_masses[0], 2)
            m_c2 = np.power(self.fs_masses[1], 2)

            m_ab2 = kin[2]
            m_ac2 = s_k_pi_wrong
            m_bc2 = kin[1]
        else:  # self.isobar_index == 23:  # ((Ks, piWrong), piRight)
            m_a2 = np.power(self.fs_masses[1], 2)
            m_b2 = np.power(self.fs_masses[2], 2)
            m_c2 = np.power(self.fs_masses[0], 2)

            m_ab2 = s_k_pi_wrong
            m_ac2 = kin[1]
            m_bc2 = kin[2]
        return m_a2, m_b2, m_c2, m_ab2, m_ac2, m_bc2


class BelleS(BelleAngDep):
    """S wave angular function."""

    def eval(self, kin):
        """Evaluate the function.

        The Breit-Wigner is always called, with the invariant mass-square of the two-particle combination, that forms
        the Isobar.
        """
        if self._check_condition_kinematic(kin):
            return complex(1., 0.)
        else:
            return complex(0., 0.)


class BelleP(BelleAngDep):
    """P wave angular function."""

    def eval(self, kin):
        """Evaluate the function.

        The Breit-Wigner is always called, with the invariant mass-square of the two-particle combination, that forms
        the Isobar.
        """
        self._check_condition_kinematic(kin)

        m_d2 = kin[0]
        s_k_pi_wrong = self.compute_s_kpi_wrong(m_d2, kin)
        m_a2, m_b2, m_c2, m_ab2, m_ac2, m_bc2 = self.compute_masses(kin, s_k_pi_wrong)

        return_value = m_ac2 - m_bc2 + (m_d2-m_c2)*(m_b2-m_a2)/m_ab2

        return complex(return_value, 0.)


class BelleD(BelleAngDep):
    """D wave angular function."""

    @staticmethod
    def _compute_return_value(m_d2, m_a2, m_b2, m_c2, m_ab2, m_ac2, m_bc2):
        return_value = m_ab2 - 2 * m_d2 - 2 * m_c2 + np.power(m_d2 - m_c2, 2)/m_ab2
        return_value *= m_ab2 - 2 * m_a2 - 2 * m_b2 + np.power(m_a2 - m_b2, 2)/m_ab2

        return_value /= -3
        return_value += np.power(m_bc2 - m_ac2 + (m_d2 - m_c2)*(m_a2 - m_b2)/m_ab2, 2)
        return return_value

    def eval(self, kin):
        """Evaluate the function.

        The Breit-Wigner is always called, with the invariant mass-square of the two-particle combination, that forms
        the Isobar.
        """
        self._check_condition_kinematic(kin)

        m_d2 = kin[0]
        s_k_pi_wrong = self.compute_s_kpi_wrong(m_d2, kin)
        m_a2, m_b2, m_c2, m_ab2, m_ac2, m_bc2 = self.compute_masses(kin, s_k_pi_wrong)

        return_value = self._compute_return_value(m_d2, m_a2, m_b2, m_c2, m_ab2, m_ac2, m_bc2)

        return complex(return_value, 0.)


class Amplitude:
    """Class implementing the amplitude as a combinaton of Breit Wigner and angular funcitons."""

    def __init__(self, mother_mass, fs_masses, identical_particles=True):

        self.identical_particles = identical_particles

        # ToDo: what is the difference between mass and mother mass?
        self.mass = mother_mass ** 2
        self.mother_mass = mother_mass
        self.s = mother_mass ** 2

        self.fs_masses = fs_masses

    def eval(self, grid):
        """Evaluate the function.

        Parameters
        ----------
        grid: ndarray
            Array containing Dalitz points as pairs, e.g. (m12, m23)

        Returns
        -------
        out: ndarray
            Array containing the evaluated output from the input Dalitz points grid.
        """

        name = 'Test'
        width = 1
        spin = 0
        rr = 0
        rd = 0

        return_value = list()
        for pair in grid:
            # pair represents the dalitz plot point to be used, eg. as m12 m13.
            bachelor_mass = self.mother_mass - (pair[0] + pair[1])

            kin = np.array([np.power(self.mother_mass, 2), pair[0] ** 2, pair[1] ** 2])

            breit_wigner_12 = BreitWigner(name, self.mass, width, spin, self.mother_mass, bachelor_mass,
                                          pair[0], pair[1], rr, rd)

            partial_wave_12 = BelleS(isobar_index=12, fs_masses=self.fs_masses)

            breit_wigner_13 = BreitWigner(name, self.mass, width, spin, self.mother_mass, bachelor_mass,
                                          pair[1], pair[0], rr, rd)

            partial_wave_13 = BelleS(isobar_index=13, fs_masses=self.fs_masses)

            if self.identical_particles:
                val = breit_wigner_12.eval(self.s) * partial_wave_12.eval(kin) \
                      + breit_wigner_13.eval(self.s) * partial_wave_13.eval(kin)
            else:
                # ToDo: what if the paticles are not identical?
                val = None

            return_value.append(val)

        return np.asarray(return_value)
