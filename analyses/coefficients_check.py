from cmath import phase as compute_phase
import matplotlib.pyplot as plt
import numpy as np

resonances_list = [ 'S:kappa', 'S:K*(1430)', 'S:const', 'P:K*(892)', 'P:K*(1410)', 'P:K*(1680)', 'D:K2*(1430)']

number_of_bins_computations = [50, 100, 200, 300, 1000]


coefficients_list_bins_50 = [-5.36340784e+02, 3.63852733e+02, 2.13197234e+02, -6.93575897e+02, -2.06161974e+02,  1.91704005e+02, -1.21226798e+02, -1.70623448e+02, 1.21149076e+02, 3.93783437e+01, -2.25384459e-01, -1.47389342e+02, 7.14309004e+01, -1.58964980e+02]

coefficients_list_bins_100 = [-558.4496132, 305.05398471, 268.14874572, -652.13551317, -226.93493899, 155.51698107, 91.2966441, 175.91360078, -113.13131799, -49.62381546, -25.81242329, 148.71431974, 82.87857374, -148.88239754]

coefficients_list_bins_200 = [-319.25314061, 264.97411517, 388.86334456, -326.8635696, 96.05808542, -74.7284627, 157.15273816,  189.63642248, -4.48689041, -4.63907273, 229.0856274,   275.86960612,  131.15062522, -108.1863837 ]

coefficients_list_bins_300 = [-289.61182997, 264.74821982, 358.45311188, -329.39457424, 83.35840405, -82.46558052, -166.48729614, -174.09935418, 2.97900924, 4.44153073, -230.00961129, -246.54968544, 120.68322997, -110.16154418]

# coefficients_list_bins_1000 = [48.90661787, 80.46677399, -174.96460306, -125.28365442,  -20.48360224, 3.02364197,  -29.91080483,  125.19118007,   96.78096618,  -98.83134796, -119.93139995,  128.63563737, -33.33647736,  -22.20312713]


class CheckCoefficientsTest:

    def __init__(self):

        fig, self.axes = plt.subplots(figsize=(6, 6))

        # self.axes.axvline(c='grey', lw=1)
        # self.axes.axhline(c='grey', lw=1)

        self.phases_list = list()


    def compute_complex_coefficients(self, coefficients):
        n = len(coefficients)

        output_list = list()
        phases_list = list()

        for i in range(int(n / 2)):
            real_part = coefficients[i * 2]
            imag_part = coefficients[i * 2 + 1]

            phase = compute_phase(real_part + 1j*imag_part)
            if phase < 0:
                phase += np.pi

            norm = np.sqrt(real_part ** 2 + imag_part ** 2)

            output_list.append((real_part + 1.0j * imag_part)/norm)
            phases_list.append(phase)

        self.phases_list.append(phases_list)

        return output_list


    def check_coefficients(self, coefficients):

        self.values = self.compute_complex_coefficients(coefficients)



    def plot(self):

        for i in range(len(self.values)):
            self.axes.scatter(self.values[i].real, self.values[i].imag, label=resonances_list[i])

        self.axes.legend()
        plt.show()

    def seperate_phases(self):
        resonance_phase = list()
        for phase_list in self.phases_list:
            resonance_phase.append(phase_list[0])
        return resonance_phase

    def plot_phases(self):
        # fig, subplots = plt.subplots(nrows=7, ncols=1, sharex='row', figsize=(12, 6))
        # for subplot in subplots:
        #     subplot.plot()

        resonance_phase = self.seperate_phases()

        plt.plot(number_of_bins_computations, resonance_phase, '*')
        plt.show()

    def plot_phases_2(self):
        for phases_list in self.phases_list:
            plt.plot(resonances_list, phases_list, '*')

        plt.xlabel('Resonance')
        plt.ylabel('Coefficient phase [radians]')

        plt.title('Coefficients phase for different number of bins')

        plt.legend(number_of_bins_computations)
        plt.show()


if __name__ == "__main__":

    check_coefficients_test = CheckCoefficientsTest()

    check_coefficients_test.check_coefficients(coefficients_list_bins_50)
    check_coefficients_test.check_coefficients(coefficients_list_bins_100)
    check_coefficients_test.check_coefficients(coefficients_list_bins_200)
    check_coefficients_test.check_coefficients(coefficients_list_bins_300)

    # check_coefficients_test.plot()
    # check_coefficients_test.plot_phases()
    check_coefficients_test.plot_phases_2()
