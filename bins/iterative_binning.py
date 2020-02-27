import numpy as np
from bins.bins_main import RectangularBin
from amplitudes.constants import m_dc, m_kc, m_pi, fs_masses
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
