#!/usr/bin/python
# chi2model.py
# Created: 2019-11-22 09:11:21.158962
# Author: Fabian Krinner

import numpy as np
import ROOT

from bins.utils import is_valid_dalitz_point, hub  # , abs2

from bins.bins_main import RectangularBin


class DalitzChi2model:
    """Class implementing a Dalitz chi2 model."""

    def __init__(self, m_mother, fs_masses, bins_list, mesh_width, functions, n_term_in_bin=10):
        """

        Parameters
        ----------
        m_mother: float
            mass of the mother particle
        fs_masses: ndarray, list[float]
            list of final state masses: [m1, m2, m3] with m12**2 on the X-axis and m13**2 on the Y-axis
        bins_list:
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
        self.number_of_bins = len(bins_list)
        self.bins_list = bins_list

        # Functions
        self.functions = functions
        self.number_of_functions = len(functions)

        # Options
        self.meshWidth = mesh_width
        self.n_term_in_bin = n_term_in_bin

        # Kinematic validity of bins
        self.number_valid_bins = None
        self.valid_bins = None  # List of kinematically allowed bins

        # Data
        self.data = None
        self.errors = None
        self.inv_errors = None

        # Integrals
        self.total_size_mesh = None
        self.number_of_mesh_bin = None
        self.integrals = None
        self.integral_indices = None
        self.norms = None

        if not self.check_contiguous(True):
            print(
                "dalitzChi2model.__init__(...): WARNING: The received bins do not contiguously cover the Dalitz plot.")

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
        for bin in self.bins_list:
            all_borders = bin.get_all_borders()
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
        @param g second function index
        """
        return f * self.number_of_functions + g

    def get_single_function_indices(self, t):
        """
        Allows one-dimensional indexing of two functions
        @param t combined function index
        """
        f = int(t / self.number_of_functions)
        g = t - self.number_of_functions * f
        return f, g

    def find_bin(self, m2_12, m2_13, break_found=True):
        """
        Find the bin index to a given point on the Dalitz Plot
        @param m2_12 invariant mass square of particles 12
        @param m2_13 invariant mass square of particles 13
        @param break_found flag if break after first found bin (If False, returns the last valid bin found)
        """
        n_found = 0
        bin_number = None
        for _bin_number, binn in enumerate(self.bins_list):
            if binn.contains(m2_12, m2_13):
                n_found += 1
                bin_number = _bin_number
                if break_found:
                    break
        return bin_number, n_found

    def check_contiguous(self, do=True):
        """
        Checks, if the received binning is contiguous
        """
        ret_val = True

        if do:
            print("Checking contigous.")
            n_grid = int(self.s / self.meshWidth)
            lin_grid = np.array([i * self.meshWidth for i in range(n_grid + 1)])
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
        return self.bins_list[n_bin].make_grid(self.meshWidth)

    def make_valid_bins(self):
        """Produce the list of valid bins i.e. the list of bins, where at least one integral point is valid."""
        self.valid_bins = np.zeros(self.number_of_bins, dtype=bool)
        print("Making valid bins.")
        for bin_number in range(self.number_of_bins):
            grid = self.make_grid(bin_number)

            valid = is_valid_dalitz_point(grid, self.s, self.fsSquare[0], self.fsSquare[1], self.fsSquare[2])

            # At least one integral point is valid
            if np.sum(valid) > 0:
                self.valid_bins[bin_number] = True

        self.number_valid_bins = np.sum(self.valid_bins)

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
            self.integral_indices[i] = np.resize(self.integral_indices[i], self.n_term_in_bin + 1)
        self.n_term_in_bin += 1

    def make_integrals(self):
        """Make the integral for all interference terms."""
        print("Making integrals.")
        self.integrals = []
        self.integral_indices = []
        self.total_size_mesh = 0
        self.number_of_mesh_bin = np.zeros(self.number_of_bins, dtype=int)
        self.norms = np.zeros(self.number_of_functions)

        for t in range(self.number_of_bins):
            print(f"Computing for bin {t}")
            if not self.valid_bins[t]:
                continue

            indices, integral, n_mesh = self.make_bin_integrals(t)

            self.total_size_mesh += n_mesh
            self.number_of_mesh_bin[t] = n_mesh
            self.integrals.append(integral)
            self.integral_indices.append(indices)
            for i, I in enumerate(indices):
                f, g = self.get_single_function_indices(I)
                if f == g:
                    self.norms[f] += integral[i]

        # Finally do the normalization
        for i, I in enumerate(self.integral_indices):
            for j, J in enumerate(I):
                f, g = self.get_single_function_indices(J)
                self.integrals[i][j] /= (self.norms[f] * self.norms[g]) ** .5

        print("Finished making integrals.")

    def make_bin_integrals(self, n_bin):
        """Make the integral for all interference terms for a single bin.

        @param n_bin bin index to calculate the integral for
        """
        if not self.valid_bins[n_bin]:
            return None

        grid = self.make_grid(n_bin)
        valid_grid = is_valid_dalitz_point(grid, self.s, self.fsSquare[0], self.fsSquare[1], self.fsSquare[2])
        n_mesh = np.sum(valid_grid)
        grid = grid[valid_grid, :]

        f_indices, lintegral = self._compute_integrals(grid)
        return f_indices, lintegral, n_mesh

    def _compute_integrals(self, grid):

        f_evals = np.array([f(grid) for f in self.functions])
        integral = np.dot(np.conjugate(f_evals), f_evals.T)
        f_indices = np.full(self.n_term_in_bin, -1, dtype=int)
        lintegral = np.zeros(self.n_term_in_bin, dtype=complex)

        count = 0
        for f in range(self.number_of_functions):
            for g in range(f, self.number_of_functions):
                if integral[f, g] == 0. + 0.j:
                    continue
                if count == self.n_term_in_bin:
                    # Do not set the last value here, since it will be overwritten in any case (won't stay empty)
                    f_indices.resize((self.n_term_in_bin + 1,))
                    lintegral.resize((self.n_term_in_bin + 1,))

                    self.increase_integral_size()

                sfi = self.get_overall_function_index(f, g)

                f_indices[count] = sfi
                lintegral[count] = integral[f, g]
                count += 1

        return f_indices, lintegral

    def get_intensity(self, prod_amps):
        """Get the model intensity for all valid bins.

        @param prod_amps production amplitudes as complex numpy array
        """
        prod_prod_amps = np.zeros(self.number_of_functions ** 2)

        self._compute_intensities(prod_amps, prod_prod_amps)
        values = self._compute_intensities_2(prod_prod_amps)
        return values

    def _compute_intensities(self, prod_amps, prod_prod_amps):
        for f in range(self.number_of_functions):
            for g in range(f, self.number_of_functions):
                sfi = self.get_overall_function_index(f, g)
                prod_prod_amps[sfi] = np.conjugate(prod_amps[f]) * prod_amps[g]
                if f == g:
                    prod_prod_amps[
                        sfi] /= 2  # Since below, we sum over 2*Re(...), the diagonal part need a factor of 1/2

    def _compute_intensities_2(self, prod_prod_amps):
        ret_val = np.zeros(self.number_valid_bins, dtype=float)

        for t in range(self.number_valid_bins):
            val = 2 * np.dot(prod_prod_amps[self.integral_indices[t]], self.integrals[t]).real
            ret_val[t] = val

        return ret_val

    def eval(self, params):
        """
        Evaluate the chi2 function
        @param params production amplitutdes as real-values numpy array [re_0, im_0, re_1, ...]
        """
        prod_amps = params[::2] + 1.j * params[1::2]
        model = self.get_intensity(prod_amps)
        return np.sum((model - self.data) ** 2 * self.inv_errors)

    def _check_number_of_data_points(self, m2_12, m2_13):
        number_data_points = len(m2_12)
        if not len(m2_13) == number_data_points:
            raise ValueError("Number of data points does not match.")
        return number_data_points

    def load_data(self, m2_12, m2_13):
        """Load the data points.

        @param m2_12 numpy array for kinematic values for m2_12 of the events
        @param m2_13 numpy array for kinematic values for m2_13 of the events
        """
        number_data_points = self._check_number_of_data_points(m2_12, m2_13)

        data = self._load_data(number_data_points, m2_12, m2_13)

        self._check_if_data_is_in_invalid_bins(data)

    def _load_data(self, number_data_points, m2_12, m2_13):
        data = np.zeros(self.number_of_bins, dtype=float)

        for data_index in range(number_data_points):
            bin_number, found_in_n_bins = self.find_bin(m2_12[data_index], m2_13[data_index], True)
            if found_in_n_bins == 0:
                error = f"dalitzChi2model.loadData(...): " \
                        f"WARNING: Could not find a valid bin for the data point: ({str(m2_12[data_index])} {str(m2_13[data_index])})."
                print(error)
                continue
            data[bin_number] += 1.

        return data

    def _check_if_data_is_in_invalid_bins(self, data):
        count_valid = 0
        self.data = np.zeros(self.number_valid_bins, dtype=float)

        for bin_index in range(self.number_of_bins):
            if not self.valid_bins[bin_index]:
                if data[bin_index] > 0:
                    error = f"dalitzChi2model.loadData(...): " \
                            f"WARNING: {str(data[bin_index])} events laying in invalid bin {str(bin_index)}: {str(self.bins_list[bin_index])}"
                    print(error)
                continue
            self.data[count_valid] = data[bin_index]
            count_valid += 1

        self.inv_errors = 1. / self.data
        self.inv_errors[np.isinf(self.inv_errors)] = 0.

    def make_theo_hist(self, params):
        """Produce a histogram for the theory.

        @params params production amplitudes to be used (real-valued numpy array: [re_0, im_0, re_1, ...])
        """
        prod_amps = params[::2] + 1.j * params[1::2]
        return self.make_hist(self.get_intensity(prod_amps), "_intens_theo")

    def make_data_hist(self):
        """Produce a histogram for the data points."""
        return self.make_hist(self.data, "_data")

    def make_hist(self, data, tag=''):
        """Produces a histogram from given data.

        @param data numpy array of the data to be plotted
        @param tag  string to name the histogram after
        """
        nx = len(self.binning_x)
        ny = len(self.binning_y)

        hist = ROOT.TH2D(f"Dalitz{tag}", f"Dalitz{tag}", nx - 1, self.binning_x, ny - 1, self.binning_y)
        delta = self.meshWidth  # offset from the actual border
        count_valid = 0
        for b in range(self.number_of_bins):
            if self.valid_bins[b]:
                if not isinstance(self.bins_list[b], RectangularBin):
                    print("WARNING: The plotting is not set up to handle non-rectangular bins.")

                borders = self.bins_list[b].get_borders()
                i_min = hist.GetXaxis().FindBin(borders[0] + delta)
                i_max = hist.GetXaxis().FindBin(borders[1] - delta)
                j_min = hist.GetYaxis().FindBin(borders[2] + delta)
                j_max = hist.GetYaxis().FindBin(borders[3] - delta)
                for i in range(i_min, i_max + 1):
                    for j in range(j_min, j_max + 1):
                        hist.SetBinContent(i, j, data[count_valid] / self.number_of_mesh_bin[b])
                count_valid += 1
        return hist
