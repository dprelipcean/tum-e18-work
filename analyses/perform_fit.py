import numpy as np


def perform_fit(c2model, fcns):
    params = np.random.uniform(-1., 1., 2 * len(fcns))
    # params = np.zeros(2 * len(fcns))

    from iminuit import Minuit

    m = Minuit.from_array_func(c2model.eval, params, error=0.5)
    print(f"Fit process started.")
    m.migrad()

    vals = np.array([m.values['x' + str(i)] for i in range(2 * len(fcns))])

    print(params)
    print(m.values)
    print(vals)

    ntfrr = ((vals[0] + 1.j * vals[1]) * (vals[2] - 1.j * vals[3])) ** 2
    print(f"{ntfrr / abs(ntfrr)} should be real (phase of pi/2)")

    return vals
