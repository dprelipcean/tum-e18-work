#!/usr/bin/python
# utils.py
# Created: 2019-11-22 09:12:34.548459
# Author: Fabian Krinner

import numpy as np


def uesen(g):

	return g[:, 0]**2


def osh(g):
	return np.zeros(g.shape[0]) + 1.1


def heb(g):
	return np.full(g.shape[0], 1.)


def hub(g):
	return g[:, 0] + g[:, 1]


def hob(g):
	return g[:, 0] * g[:, 1]


def half_budalf(g):
	ret_val = g[:, 0] + g[:, 1]
	ret_val[g[:, 0] < 1.] = 0.
	return ret_val


def abs2(x):
	return x.real**2 + x.imag**2


def q2(m2, m12, m22):
	return (m2 ** 2 + m12 ** 2 + m22 ** 2 - 2. * (m2 * m12 + m2 * m22 + m12 * m22)) / 4 / m2


def is_valid_dalitz_point(dalitz_points, s, m12, m22, m32):
	"""Check whether the given dalitz point is valid.

	Parameters
	----------
	dalitz_points: ndarray
	s: float
		Centre of mass energy.
	m12: float
	m22: float
	m32: float

	Returns
	-------

	"""

	s23 = s + m12 + m22 + m32 - dalitz_points[:, 0] - dalitz_points[:, 1]

	ret_val = s23 >= (m22**.5+m32**.5)**2
	ret_val &= s23 <= (s**.5-m12**.5)**2
	p12 = q2(s, s23, m12)
	p22 = q2(s, dalitz_points[:, 1], m22)
	p32 = q2(s, dalitz_points[:, 0], m32)
	q12 = q2(dalitz_points[:, 0], m12, m22)
	q13 = q2(dalitz_points[:, 1], m12, m32)
	q23 = q2(s23, m22, m32)

	ret_val &= p12 >= 0.
	ret_val &= p22 >= 0.
	ret_val &= p32 >= 0.
	ret_val &= q12 >= 0.
	ret_val &= q13 >= 0.
	ret_val &= q23 >= 0.

	p1 = abs(p12)**.5  # can take the abs() here, since p < 0. is already set to false in retVal
	p2 = abs(p22)**.5
	p3 = abs(p32)**.5

	ret_val &= p1 + p2 > p3
	ret_val &= p1 + p3 > p2
	ret_val &= p2 + p3 > p1

	return ret_val.flatten()


def main():
	import ROOT
	m_dc = 1.86958
	m_pi = 0.13957
	m_kc = 0.493677
	n_bins = 1000
	hist = ROOT.TH2D("f", "f", n_bins, 0., m_dc**2, n_bins, 0., m_dc**2)
	bin_centers = []
	for i in range(n_bins):
		for j in range(n_bins):
			bin_centers.append([hist.GetXaxis().GetBinCenter(i+1), hist.GetYaxis().GetBinCenter(j+1)])
	bin_centers = np.array(bin_centers)

	valid = is_valid_dalitz_point(bin_centers, m_dc ** 2, m_kc ** 2, m_pi ** 2, m_pi ** 2)
	count = 0
	for v in valid:
		if v:
			count += 1
	print(f'sdfsfa{count}')
	centers_to_set = bin_centers[valid, :]
	print(centers_to_set.shape)
	for cts in centers_to_set:
		i = hist.GetXaxis().FindBin(cts[0])
		j = hist.GetYaxis().FindBin(cts[1])
		hist.SetBinContent(i, j, 1.)
	hist.Draw("col")
	input()


if __name__ == "__main__":
	main()