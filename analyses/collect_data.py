import numpy as np


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
    m2s = np.random.uniform((m_pi + m_kc) ** 2, (m_dc - m_pi) ** 2, 2 * n_data)
    m2s.shape = (n_data, 2)
    return m2s


def read_data_monte_carlo(data_range=None):
    """Read the monte carlo generated data file."""
    with open('../data/CP_MC_data_SPD.CP_MC') as f:
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

    data_size = len(data)
    if data_range:
        data_array = np.asarray(data[0:data_range])
        data_array.shape = (int(data_range / 2), 2)

    else:
        data_array = np.asarray(data)
        data_array.shape = (int(data_size / 2), 2)

    return data_array, data_size
