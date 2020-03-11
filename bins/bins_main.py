#!/usr/bin/python
# bins_main.py
# Created: 2020-01-03 09:10:40.360767
# Author: Fabian Krinner
# Edit by: Daniel Prelipcean

"""Script that defines bins functions for the Dalitz plot."""

import numpy as np


class RectangularBin:
    """Class implementing rectangular bins."""

    def __init__(self, x_min, x_max, y_min, y_max):
        self.borders = [x_min, x_max, y_min, y_max]
        self.value = 0

        self.position = ((x_min + x_max)/2, (y_min + y_max)/2)

    def get_all_borders(self):
        """Return the borders of the bins as a list."""
        return [self.borders]

    def get_borders(self):
        """Return the borders of the bins."""
        return self.borders

    def increment_value(self):
        """Increment the bin value by one."""
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
        return self.borders[0] <= x <= self.borders[1] and self.borders[2] <= y <= self.borders[3]

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
    """Class that creates bins."""

    def __init__(self, m_dc, n_bins, bin_class=RectangularBin, isotropic_binning=True):
        self.bin_class = bin_class

        self.m_dc = m_dc
        self.n_bins = n_bins

        self.resonance_regions = [[0.8, 1.2], [1.2, 0.8], [1.6, 1.6]]
        self.central_resonance_region = 1.75
        self.division_factor = n_bins/10

        self.binning_x = self._create_binning_1d()
        self.dx = self._get_bin_increment(self.binning_x)

        self.binning_y = self._create_binning_1d()
        self.dy = self._get_bin_increment(self.binning_y)

        self.isotropic_binning = isotropic_binning

    def _create_binning_1d(self):
        return np.linspace(0., self.m_dc**2, self.n_bins)

    @staticmethod
    def _get_bin_increment(binning_1d):
        dx = binning_1d[1] - binning_1d[0]
        return dx

    def _compute_distance_from_resonances(self, position, ds, position2, ds2):
        distance_from_resonance = self.m_dc ** 2/ ds
        # print(distance_from_resonance)
        individual_contributions = [0, 0]

        for resonance in self.resonance_regions:
            distance_from_resonance_new, individual_contributions_new = self._measure_function(position, resonance[0], position2, resonance[1])

            if distance_from_resonance_new < distance_from_resonance:
                distance_from_resonance = distance_from_resonance_new
                individual_contributions = individual_contributions_new
        # print(f'{distance_from_resonance} {individual_contributions}')

        return individual_contributions

    def _compute_distance_from_resonance_center(self, _bin, ds):
        distance_from_resonance_center = abs(_bin - self.central_resonance_region) / ds
        return distance_from_resonance_center

    def _compute_distance_increment(self, distance_from_resonance):
        if self.isotropic_binning:
            distance_increment = 1
        else:
            distance_increment = round(distance_from_resonance / self.division_factor)

            if distance_increment == 0:
                distance_increment = 1

        return distance_increment

    def _compute_end_index(self, index, distance_increment):
        end_index = index + distance_increment
        if end_index >= self.n_bins:
            end_index = self.n_bins - 1

        return end_index

    def _create_individual_bin(self, i, end_index_i, j, end_index_j):
        return self.bin_class(self.binning_x[i], self.binning_x[int(end_index_i)],
                              self.binning_y[j], self.binning_y[int(end_index_j)])

    def _measure_function(self, x, x0, y, y0):
        x_contribution = abs(x - x0) / self.dx
        y_contribution = abs(y - y0) / self.dy
        contribution = np.sqrt(x_contribution ** 2 + y_contribution ** 2)
        return contribution, (x_contribution, y_contribution)

    def _compute_increments(self, index, binning, ds, index2, binning2, ds2):
        distance_from_resonance, distance_from_resonance2 = self._compute_distance_from_resonances(binning[index], ds, binning2[index2], ds2)

        # distance_from_resonance_center = self._compute_distance_from_resonance_center(binning[index], ds)
        # distance_from_resonance_center2 = self._compute_distance_from_resonance_center(binning2[index2], ds2)
        # if distance_from_resonance_center < distance_from_resonance \
        #         and distance_from_resonance_center2 < distance_from_resonance2:
        #     distance_from_resonance = min(distance_from_resonance, distance_from_resonance_center)
        #     distance_from_resonance2 = min(distance_from_resonance2, distance_from_resonance_center2)

        distance_increment = self._compute_distance_increment(distance_from_resonance)
        distance_increment2 = self._compute_distance_increment(distance_from_resonance2)

        return distance_increment, distance_increment2

    def create_bins(self):
        """Create the bins for tor the Dalitz plot."""

        bin_list = []

        i = 0
        while i < (self.n_bins - 1):
            end_index_i = i
            j = 0
            while j < (self.n_bins - 1):
                distance_increment_i, distance_increment_j \
                    = self._compute_increments(i, self.binning_x, self.dx, j, self.binning_y, self.dy)

                end_index_i = self._compute_end_index(i, distance_increment_i)
                end_index_j = self._compute_end_index(j, distance_increment_j)

                # print(f'{i} {j} {distance_increment_i} {distance_increment_j} {end_index_i} {end_index_j}')

                bin_list.append(self._create_individual_bin(i, end_index_i, j, end_index_j))

                j = int(end_index_j)
            i = int(end_index_i)

        # for bin in bin_list:
        #     print(bin.get_borders())
        return bin_list, (self.binning_x, self.binning_y)


def create_bins(m_dc, n_bins):
    """Create explicit bins."""
    bins = Bins(m_dc, n_bins)
    return bins.create_bins()
