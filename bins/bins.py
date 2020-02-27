#!/usr/bin/python
# bins.py
# Created: 2020-01-03 09:10:40.360767
# Author: Fabian Krinner

"""Script that defines bins functions for the Dalitz plot."""

import numpy as np


class RectangularBin:
    """Class implementing rectangular bins."""

    def __init__(self, x_min, x_max, y_min, y_max):
        self.borders = [x_min, x_max, y_min, y_max]
        self.value = 0

    def get_all_borders(self):
        """Return the borders of the bins as a list."""
        return [self.borders]

    def get_borders(self):
        """Return the borders of the bins."""
        return self.borders

    def increment_value(self):
        self.value += 1

    def contains(self, x, y):
        """Check whether a point is contained within the bins' borders.

        Parameters
        ----------
        x: int, float
        y: int, float
            Points in the plane.

        Returns
        -------
        out: bool
            Flag indicating whether the point is withing the bins' borders or not.
        """
        return self.borders[0] < x <= self.borders[1] and self.borders[2] < y <= self.borders[3]

    def make_grid(self, mesh_width):
        """

        Parameters
        ----------
        mesh_width: int, float
            The width for the rectangular bins.

        Returns
        -------
        grid: np.array
            The bins as an array.
        """

        i_min = int(self.borders[0] / mesh_width) + 1
        i_max = int(self.borders[1] / mesh_width) + 1

        j_min = int(self.borders[2] / mesh_width) + 1
        j_max = int(self.borders[3] / mesh_width) + 1

        p_x = np.array([i * mesh_width for i in range(i_min, i_max)])
        p_y = np.array([j * mesh_width for j in range(j_min, j_max)])

        grid = np.array(np.meshgrid(p_x, p_y)).T.reshape(-1, 2)
        return grid


class Bins:

    def __init__(self, m_dc, n_bins):
        self.m_dc = m_dc
        self.n_bins = n_bins

        self.resonance_region = 0.8
        self.central_resonance_region = 1.75
        self.division_factor = n_bins/20

    def _create_binning_1d(self):
        return np.linspace(0., self.m_dc**2, self.n_bins)

    def _get_bin_increment(self, binning_1d):
        dx = binning_1d[1] - binning_1d[0]
        return dx

    def create_bins(self, m_dc, n_bins=20):
        """Create the bins for tor the Dalitz plot."""
        binning_x = self._create_binning_1d()
        dx = self._get_bin_increment(binning_x)

        binning_y = self._create_binning_1d()
        dy = self._get_bin_increment(binning_y)

        bin_list = []

        i = 0
        while i < (self.n_bins - 1):
            distance_from_resonance_x = (abs(binning_x[i] - self.resonance_region) / dx)
            increment_distance_i = round(distance_from_resonance_x / self.division_factor)

            j = 0
            while j < (self.n_bins - 1):
                distance_from_resonance_y = (abs(binning_y[j] - self.resonance_region) / dy)
                increment_distance_j = round(distance_from_resonance_y / self.division_factor)

                if increment_distance_i == 0:
                    increment_distance_i = 1
                if increment_distance_j == 0:
                    increment_distance_j = 1

                end_index_i = i + increment_distance_i
                if end_index_i >= n_bins:
                    end_index_i = n_bins - 1
                end_index_j = j + increment_distance_j
                if end_index_j >= n_bins:
                    end_index_j = n_bins - 1

                print(f'{i} {j} {increment_distance_i} {increment_distance_j} {end_index_i} {end_index_j}')
                bin_list.append(RectangularBin(binning_x[i], binning_x[int(end_index_i)],
                                               binning_y[j], binning_y[int(end_index_j)]))

                j = int(end_index_j)
            i = int(end_index_i)

        for bin in bin_list:
            print(bin.get_borders())
        return bin_list, (binning_x, binning_y)


def create_bins(m_dc, n_bins):
    bins = Bins(m_dc, n_bins)
    return bins.create_bins(m_dc, n_bins)
