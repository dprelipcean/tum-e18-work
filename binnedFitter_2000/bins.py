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

	def get_all_borders(self):
		"""Return the borders of the bins as a list."""
		return [self.borders]

	def get_borders(self):
		"""Return the borders of the bins."""
		return self.borders

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
