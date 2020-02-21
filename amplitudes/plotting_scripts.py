import matplotlib.pyplot as plt


def plot(x_points, points_real, points_imag, points_abs):
    plt.plot(x_points, points_real)
    plt.plot(x_points, points_imag)
    plt.plot(x_points, points_abs)

    plt.legend(['Real part', 'Imaginary part', 'Magnitude'])

    plt.show()
