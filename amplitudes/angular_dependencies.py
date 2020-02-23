# -*- coding: utf-8 -*-
#
# Copyright (C) 2019, 2020 TUM E18 Research Group.

"""Classes implementing the angular dependencies for the S, P and D wave functions."""

import numpy as np


class KinematicSignature:
    """Class for checking the kinematic signature."""
    nKin = 3


class AngularDependence:
    """Parent class for angular dependence functions."""

    valid_isobar_indexes = [12, 13, 23]

    def __init__(self):
        pass

    def _check_condition_isobarindex(self, isobar_index):
        """Check that the isobar index is among the valid values.

        Parameters
        ----------
        isobar_index: int
        """
        if isobar_index not in self.valid_isobar_indexes:
            raise Exception(f"None of the three possible isobar masses match (12, 13 and 23) the given value: "
                            f"{isobar_index}")

    @staticmethod
    def _check_condition_fsmasses(fs_masses):
        """Check that the fs_masses respect some conditions.

        Parameters
        ----------
        fs_masses: ndarray, tuple[float}
        """
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
        fs_masses: ndarray, list
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
