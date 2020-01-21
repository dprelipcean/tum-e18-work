#!/usr/bin/python
# chi2model.py
# Created: 2019-11-22 09:11:21.158962
# Author: Fabian Krinner

import numpy as np
import ROOT

from utils import is_valid_dalitz_point, hub  # , abs2

import bins
from ppppppppp.massShape import Amplitude


class DalitzChi2model:
	"""Class implementing a Dalitz chi2 model."""
	def __init__(self, m_mother, fs_masses, _bins, mesh_width, functions, n_term_in_bin=10):
		"""

		Parameters
		----------
		m_mother: float
			mass of the mother particle
		fs_masses: ndarray, list[float]
			list of final state masses: [m1, m2, m3] with m12**2 on the X-axis and m13**2 on the Y-axis
		_bins:
			list of bins to use. Format [(m12_low, m12_high, m13_low, m13_high), ... ]
		mesh_width: float
			distance of two points in one direction on the integral mesh
		functions:
			list of functions ( = partialWaves)
		n_term_in_bin: int, optional
			(estimated) number of non-zero terms in every bin ( = nonZeroAmpls**2)
			Defaults to 10.
		"""

		# Physics
		self.mMother = m_mother
		self.s = m_mother ** 2
		self.fsMasses = fs_masses
		self.fsSquare = [m ** 2 for m in fs_masses]

		# Binning
		self.nBins = len(_bins)
		self.bins = _bins

		# Functions
		self.functions = functions
		self.nFunc = len(functions)

		# Options
		self.meshWidth = mesh_width
		self.n_term_in_bin = n_term_in_bin

		# Kinematic validity of bins
		self.nValid = None
		self.validBins = None  # List of kinematically allowed bins

		# Data
		self.data = None
		self.errors = None
		self.inv_errors = None

		# Integrals
		self.totalNmesh = None
		self.nMeshBin = None
		self.integrals = None
		self.integralIndices = None
		self.norms = None

		if not self.check_contiguous(False):
			print("dalitzChi2model.__init__(...): WARNING: The received bins do not contiguously cover the Dalitz plot.")

		# Plotting
		self.binning_x = None
		self.binning_y = None
		self.make_maximal_binning()

	def make_maximal_binning(self):
		"""
		Creates all bin borders 
		"""
		print("Making maximal binning.")
		binning_x = []
		binning_y = []
		for b in self.bins:
			all_borders = b.get_all_borders()
			for borders in all_borders:
				if not borders[0] in binning_x:
					binning_x.append(borders[0])
				if not borders[1] in binning_x:
					binning_x.append(borders[1])
				if not borders[2] in binning_y:
					binning_y.append(borders[2])
				if not borders[3] in binning_y:
					binning_y.append(borders[3])
		binning_x.sort()
		binning_y.sort()
		self.binning_x = np.array(binning_x, dtype=np.float64)
		self.binning_y = np.array(binning_y, dtype=np.float64)
		print("Maximal binning finished.")

	def get_overall_function_index(self, f, g):
		"""
		Allows one-dimensional indexing of two functions
		@param f first function index
		@param g second funtcion index
		"""
		return f * self.nFunc + g

	def get_single_function_indices(self, t):
		"""
		Allows one-dimensional indexing of two functions
		@param t combined function index
		"""
		f = int(t/self.nFunc)
		g = t - self.nFunc * f
		return f, g

	def find_bin(self, m2_12, m2_13, break_found=True):
		"""
		Find the bin index to a given point on the Dalitz Plot
		@param m2_12 invariant mass square of particles 12
		@param m2_13 invariant mass square of particles 13
		@param break_found flag if break after first found bin (If False, returns the last valid bin found)
		"""
		n_found = 0
		n_bin = None
		for b, binn in enumerate(self.bins):
			if binn.contains(m2_12, m2_13):
				n_found += 1
				n_bin = b
				if break_found:
					break
		return n_bin, n_found

	def check_contiguous(self, do=True):
		"""
		Checks, if the received binning is contiguous
		"""
		ret_val = True

		if do:
			print("Checking contigous.")
			n_grid = int(self.s/self.meshWidth)
			lin_grid = np.array([i * self.meshWidth for i in range(n_grid+1)])
			grid = np.array(np.meshgrid(lin_grid, lin_grid)).T.reshape(-1, 2)
			grid = grid[is_valid_dalitz_point(grid, self.s, self.fsSquare[0], self.fsSquare[1], self.fsSquare[2])]
			for i in range(grid.shape[0]):
				t, n = self.find_bin(grid[i, 0], grid[i, 1], break_found=False)
				if n == 0:
					error = f'dalitzChi2model.checkContiguous(): ' \
							f'WARNING: No valid bin found for valid grid point: ({str(grid[i, 0])} {str(grid[i, 1])}).'
					print(error)
					ret_val = False
				elif n > 1:
					error = f'dalitzChi2model.checkContiguous(): ' \
							f'WARNING: {str(n)} valid bins found for valid grid point: ({str(grid[i, 0])} {str(grid[i, 1])}).'
					print(error)
					ret_val = False

			print("Finished hecking contigous.")
		else:
			print("Skipped checking contigous.")
		return ret_val
		
	def make_grid(self, n_bin):
		"""
		Makes the integral meshgrid for a 1-D bin index
		@param n_bin bin index to make the grid for
		"""
		return self.bins[n_bin].make_grid(self.meshWidth)

		# binBorders = self.bins[nBin]
		#
		# iMin = int(binBorders[0]/self.meshWidth) + 1
		# iMax = int(binBorders[1]/self.meshWidth) + 1
		#
		# jMin = int(binBorders[2]/self.meshWidth) + 1
		# jMax = int(binBorders[3]/self.meshWidth) + 1
		#
		# pX = np.array([i*self.meshWidth for i in range(iMin, iMax)])
		# pY = np.array([j*self.meshWidth for j in range(jMin, jMax)])
		#
		# grid = np.array(np.meshgrid(pX,pY)).T.reshape(-1,2)
		# return grid

	def make_valid_bins(self):
		"""
		Produces the list of valid bins i.e. the list of bins, where at least one integral point is valid
		"""
		self.validBins = np.zeros(self.nBins, dtype=bool)
		for b in range(self.nBins):
			grid = self.make_grid(b)
			valid = is_valid_dalitz_point(grid, self.s, self.fsSquare[0], self.fsSquare[1], self.fsSquare[2])
			if np.sum(valid) > 0:
				self.validBins[b] = True
		for i, v in enumerate(self.validBins):
			if v:
				print(i)

		self.nValid = np.sum(self.validBins)

	def increase_integral_size(self):
		"""
		Increases the storage size of the integral values for each bin by one
		"""
		error = f"dalitzChi2model.increaseIntegralSize(): " \
			f"WARNING: The estimated number of non-zero amplitudes ({str(self.n_term_in_bin)}) is too small. " \
			f"Increase by one."
		print(error)
		for i in range(len(self.integrals)):
			self.integrals[i] = np.resize(self.integrals[i], self.n_term_in_bin + 1)
			self.integralIndices[i] = np.resize(self.integralIndices[i], self.n_term_in_bin + 1)
		self.n_term_in_bin += 1

	def make_integrals(self):
		"""
		Make the integral for all interference terms
		"""
		print("Making integrals.")
		self.integrals = []
		self.integralIndices = []
		self.totalNmesh = 0
		self.nMeshBin = np.zeros(self.nBins, dtype=int)
		self.norms = np.zeros(self.nFunc)
		for t in range(self.nBins):
			if not self.validBins[t]:
				continue
			indices, integral, n_mesh = self.make_bin_integrals(t)
			self.totalNmesh += n_mesh
			self.nMeshBin[t] = n_mesh
			self.integrals.append(integral)
			self.integralIndices.append(indices)
			for i, I in enumerate(indices):
				f, g = self.get_single_function_indices(I)
				if f == g:
					self.norms[f] += integral[i]

		# Finally do the normalization
		for i, I in enumerate(self.integralIndices):
			for j, J in enumerate(I):
				f, g = self.get_single_function_indices(J)
				self.integrals[i][j] /= (self.norms[f] * self.norms[g])**.5
		print("Finished making integrals.")

	def make_bin_integrals(self, n_bin):
		"""
		Make the integral for all interference terms for a single bin
		@param n_bin bin index to calculate the integral for
		"""
		if not self.validBins[n_bin]:
			return None
		grid = self.make_grid(n_bin)
		# print(f'MyGrid is: {grid}')
		valid_grid = is_valid_dalitz_point(grid, self.s, self.fsSquare[0], self.fsSquare[1], self.fsSquare[2])
		n_mesh = np.sum(valid_grid)
		grid = grid[valid_grid, :]
		f_evals = np.array([f(grid) for f in self.functions])
		integral = np.dot(np.conjugate(f_evals), f_evals.T)
		f_indices = np.full(self.n_term_in_bin, -1, dtype=int)
		lintegral = np.zeros(self.n_term_in_bin, dtype=complex)
		count = 0
		for f in range(self.nFunc):
			for g in range(f, self.nFunc):
				if integral[f, g] == 0.+0.j:
					continue
				if count == self.n_term_in_bin:
					# Do not set the last value here, since it will be overwritten in any case (won't stay empty)
					f_indices.resize((self.n_term_in_bin + 1, ))
					lintegral.resize((self.n_term_in_bin + 1, ))
					self.increase_integral_size()

				sfi = self.get_overall_function_index(f, g)
				f_indices[count] = sfi
				lintegral[count] = integral[f, g]
				count += 1
		return f_indices, lintegral, n_mesh

	def get_intensity(self, prod_amps):
		"""
		Get the model intensity for all valid bins.
		@param prod_amps production amplotudes as complex numpy array
		"""
		ret_val = np.zeros(self.nValid, dtype=float)
		prod_prod_amps = np.zeros(self.nFunc**2)
		for f in range(self.nFunc):
			for g in range(f, self.nFunc):
				sfi = self.get_overall_function_index(f, g)
				prod_prod_amps[sfi] = np.conjugate(prod_amps[f]) * prod_amps[g]
				if f == g:
					prod_prod_amps[sfi] /= 2  # Since below, we sum over 2*Re(...), the diagonal part need a factor of 1/2
		for t in range(self.nValid):
			val = 2*np.dot(prod_prod_amps[self.integralIndices[t]], self.integrals[t]).real
			ret_val[t] = val
		return ret_val

	def eval(self, params):
		"""
		Evaluate the chi2 function
		@param params production amplutdes as real-values numpy array [re_0, im_0, re_1, ...]
		"""
		prod_amps = params[::2] + 1.j * params[1::2]
		model = self.get_intensity(prod_amps)
		return np.sum((model - self.data) ** 2 * self.inv_errors)

	def load_data(self, m2_12, m2_13):
		"""
		Load the data points
		@param m2_12 numpy array for kinematic values for m2_12 of the events
		@param m2_13 numpy array for kinematic values for m2_13 of the events
		"""
		n_dat = len(m2_12)
		if not len(m2_13) == n_dat:
			raise ValueError("Number of data points does not match.")
		data = np.zeros(self.nBins, dtype=float)
		for d in range(n_dat):
			t, n = self.find_bin(m2_12[d], m2_13[d], True)
			if n == 0:
				error = f"dalitzChi2model.loadData(...): " \
						f"WARNING: Could not find a valid bin for the data point: ({str(m2_12[d])} {str(m2_13[d])})."
				print(error)
				continue
			data[t] += 1.
		count_valid = 0
		self.data = np.zeros(self.nValid, dtype=float)
		for t in range(self.nBins):
			if not self.validBins[t]:
				if data[t] > 0:
					error = f"dalitzChi2model.loadData(...): " \
							f"WARNING: {str(data[t])} events laying in invalid bin {str(t)}: {str(self.bins[t])}"
					print(error)
				continue
			self.data[count_valid] = data[t]
			count_valid += 1
		self.inv_errors = 1. / self.data ** .5
		self.inv_errors[np.isinf(self.inv_errors)] = 0.

	def make_theo_hist(self, params):
		"""
		Produces a histogram for the theory
		@params params production amplitudes to be used (real-valued numpy array: [re_0, im_0, re_1, ...])
		"""
		prod_amps = params[::2] + 1.j * params[1::2]
		return self.make_hist(self.get_intensity(prod_amps), "_intens_theo")

	def make_data_hist(self):
		"""
		Produces a histogram for the data points
		"""
		return self.make_hist(self.data, "_data")

	def make_hist(self, data, tag=''):
		"""
		Produces a histogram from given data
		@param data numpy array of the data to be plotted
		@param tag  string to name the histogram after
		"""
		nx = len(self.binning_x)
		ny = len(self.binning_y)

		hist = ROOT.TH2D(f"Dalitz{tag}", f"Dalitz{tag}", nx - 1, self.binning_x, ny - 1, self.binning_y)
		delta = self.meshWidth  # offset from the actual border
		count_valid = 0
		for b in range(self.nBins):
			if self.validBins[b]:
				if not isinstance(self.bins[b], bins.RectangularBin):
					print("WARNING: The plotting is not set up to handle non-rectangular bins.")

				borders = self.bins[b].get_borders()
				i_min = hist.GetXaxis().FindBin(borders[0] + delta)
				i_max = hist.GetXaxis().FindBin(borders[1] - delta)
				j_min = hist.GetYaxis().FindBin(borders[2] + delta)
				j_max = hist.GetYaxis().FindBin(borders[3] - delta)
				for i in range(i_min, i_max+1):
					for j in range(j_min, j_max+1):
						hist.SetBinContent(i, j, data[count_valid]/self.nMeshBin[b])
				count_valid += 1
		return hist


def generate_random_data(m_pi, m_kc, m_dc, n_data):
	"""Generate random data for testing purposes.

	Parameters
	----------
	m_pi
	m_kc
	m_dc
	n_data: int
		Number of data points to generate.

	Returns
	-------
	m2s: ndarray
		N x 2 array of points.
	"""
	m2s = np.random.uniform((m_pi+m_kc)**2, (m_dc - m_pi)**2, 2*n_data)
	m2s.shape = (n_data, 2)
	return m2s


def read_data_monte_carlo():
	"""Read the monte carlo generated data file."""
	with open('../CP_MC_data_SPD.CP_MC') as f:
		data = list()
		for row in f:
			mother, bachelor, isobar = row.split()
			isobar = isobar[:-1]

			try:
				data.append(float(bachelor))
				data.append(float(isobar))
			except ValueError:
				data.pop(-1)
				continue

	a = np.asarray(data)
	a.shape = (int(len(data)/2), 2)
	return a


def create_bins(m_dc):
	"""Create the bins for tor the Dalitz plot."""
	n_bins = 20

	binning_x = np.linspace(0., m_dc**2, n_bins)
	binning_y = np.linspace(0., m_dc**2, n_bins)

	bin_list = []
	for i in range(len(binning_x)-1):
		for j in range(len(binning_y)-1):
			bin_list.append(bins.RectangularBin(binning_x[i], binning_x[i + 1], binning_y[j], binning_y[j + 1]))
	return bin_list


def main():
	m_dc = 1.86958
	m_pi = 0.13957
	m_kc = 0.493677

	fs_masses = [m_kc, m_pi, m_pi]

	bin_list = create_bins(m_dc)
	
	# fcns = [osh, heb, hub, hob, uesen, halfBudalf]
	amplitude = Amplitude(mother_mass=m_dc, fs_masses=fs_masses)
	fcns = [amplitude.eval]
	print("Amplitude created and functions assigned.")

	c2model = DalitzChi2model(m_dc, fs_masses, bin_list, 0.001, fcns, 10)
	print("Chi2 model created.")

	c2model.make_valid_bins()
	c2model.make_integrals()

	# m2s = generate_random_data(m_pi, m_kc, m_dc,< n_data=10000)
	m2s = read_data_monte_carlo()

	valid = is_valid_dalitz_point(m2s, m_dc ** 2, fs_masses[0] ** 2, fs_masses[1] ** 2, fs_masses[2] ** 2)
	print(f"nGOOD ={np.sum(valid)}")

	m2s = m2s[valid, :]
	weights = 1. + hub(m2s)**2

	mxx = np.max(weights)

	weights -= np.random.random(weights.shape) * mxx

	m2s = m2s[weights > 0., :]

	print(f"nDeWeight ={np.sum(weights > 0)}")

	c2model.load_data(m2s[:, 0], m2s[:, 1])

	# params = np.random.uniform(-1.,1.,2*len(fcns))
	params = np.zeros(2*len(fcns))		

	from iminuit import Minuit

	m = Minuit.from_array_func(c2model.eval, params, error=0.5)
	m.migrad()

	vals = np.array([m.values['x' + str(i)] for i in range(2*len(fcns))])

	h1 = c2model.make_theo_hist(vals)
	h1.Draw('col')
	input()
	h2 = c2model.make_data_hist()
	h2.Draw('col')
	input()
	print(vals)

	ntfrr = ((vals[0] + 1.j*vals[1]) * (vals[2] - 1.j * vals[3])) ** 2
	print(f"{ntfrr/abs(ntfrr)} should be real (phase of pi/2)")


if __name__ == "__main__":
	main()
