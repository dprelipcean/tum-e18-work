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


if __name__ == '__main__':
    spline = SplineFunction([0, 1])
    spline.plot()
