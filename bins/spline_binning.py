import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate


def step_function(x, x_min=0, x_max=1):
    if isinstance(x, np.ndarray):
        return_list = list()
        for value in x:
            return_list.append(step_function(value))
        return return_list

    else:
        if x_min <= x <= x_max:
            return_value = 1
        else:
            return_value = 0
        return return_value


class SplineFunction:

    def __init__(self, bin_list):
        self.bin_list = bin_list

        self.n = len(self.bin_list)

    def _decide_y_points(self, spline_degree):
        if spline_degree == 0:
            y_points = np.array([0, 1])
        elif spline_degree == 1:
            y_points = np.array([0., 1., 0.])
        elif spline_degree == 2:
            y_points = np.array([0., 0.9, 0.9, 0.])
        elif spline_degree == 3:
            y_points = np.array([0., 0.75, 1, 0.75, 0.])
        elif spline_degree == 4:
            y_points = np.array([0., 0.75, 1, 1, 0.75, 0.])
        elif spline_degree == 5:
            y_points = np.array([0., 0.25, 0.75, 1, 0.75, 0.25, 0.])

        return y_points

    def _define_spline(self, spline_degree):
        x_points = np.linspace(0, 1, spline_degree+2)

        y_points = self._decide_y_points(spline_degree)
        print(f'{x_points} {y_points}')
        self.spline_coefficients = interpolate.splrep(x_points, y_points, k=spline_degree)

    def call(self):
        pass

    def plot(self):
        x_points = np.linspace(-0.1, 1.1, 100)

        fig, axes = plt.subplots(1, 6)

        axes[0].plot(x_points, step_function(x_points))

        for spline_degree in [1, 2, 3, 4, 5]:
            self._define_spline(spline_degree)
            spline_values = interpolate.splev(x_points, self.spline_coefficients)

            axes[spline_degree].plot(x_points, spline_values)

        plt.show()

def draft_main():
    spline = SplineFunction([0, 1])
    spline.plot()

class Spline:

    def __init__(self, s_min, s_max, y_val, deg):
        self.s_min = s_min
        self.s_max = s_max

        self.y_val = y_val

        self.deg = deg

    def eval(self, x, nDiff=0):
        x = (x - self.s_min) / (self.s_max - self.s_min)

        if x <= 0. or x > 1.:
            return 0.

        if nDiff == 0:
            if self.deg == 0:
                return_value = 1
            elif self.deg == 1:
                if x < 0.5:
                    return_value = 2 * x
                else:
                    return_value = 2 - 2 * x
            elif self.deg == 2:
                return_value = x ** 2 - x + 1. / 6.
            return_value *= self.y_val
        elif nDiff == 1:
            if self.deg == 1:
                if x < 0.5:
                    return_value = 2
                else:
                    return_value = - 2

        return return_value


def compute_smoothness():
    pass


def spline_bining():
    x_binning_values = [0, 1, 2, 3, 4, 5, 6]
    y_values = [0, 1, 3, 2, 4, 2, 0]

    splines = list()

    for i in range(len(x_binning_values) - 1):
        s_min = x_binning_values[i]
        s_max = x_binning_values[i + 1]

        spline = Spline(s_min, s_max, y_values[i], deg=2)

        splines.append(spline)

    x_points = np.linspace(x_binning_values[0], x_binning_values[-1], 100)
    y_points = list()

    for point in x_points:
        val = 0
        for spline in splines:
            val += spline.eval(point)
        y_points.append(val)

    print(smoothness_difference(x_binning_values, splines))

    plt.plot(x_points, y_points)
    plt.show()


def smoothness_difference(x_points, splines):
    diff = 0
    for i in range(len(x_points)-1):
        for j in range(len(splines)-1):
            diff += abs(splines[j].eval(x_points[i + 1]) - splines[j+1].eval(x_points[i+1]))
    return diff


if __name__ == '__main__':
    # draft_main()
    spline_bining()
