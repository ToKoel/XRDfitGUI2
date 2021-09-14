import numpy as np
import numexpr as ne
from lmfit import minimize, Parameters, fit_report


class Data:
    def __init__(self, file, wavelength, bkg=None):
        self.file = file
        try:
            self.x, self.y, self.ey = np.loadtxt(file, unpack=True)
            if bkg:
                self.bkg_x, self.bkg_y, self.bkg_ey = np.loadtxt(bkg, unpack=True)
            else:
                self.bkg_x, self.bkg_y, self.bkg_ey = 0, 0, 0
        except ValueError:
            self.x, self.y = np.loadtxt(file, unpack=True, skiprows=4)
            if bkg:
                self.bkg_x, self.bkg_y = np.loadtxt(bkg, unpack=True, skiprows=4)
            else:
                self.bkg_x, self.bkg_y = 0, 0

        self.Q = 4 * np.pi / wavelength * np.sin(self.x / 2 * np.pi / 180)
        self.intensity = (self.y - self.bkg_y) / 50


class Fit:
    def __init__(self, q, I, param_dict=None, name=""):
        self.q = q
        self.I = I
        self.name = name
        self.params = Parameters()
        self.bkg_poly = None
        self.result_params = self.params
        self.fit_curve = np.zeros_like(q)

    def update_parameters(self, param_dict):
        for key,value in param_dict.items():
            if key == 'bkgt':
                self.params.add(key, value=value[0], min=-0.2, vary=value[1])
            else:
                self.params.add(key, value=value[0], min=0, vary=value[1])
        self.params.pretty_print()

    def fit_function(self, var_params):
        positions = np.array([var_params['P1'], var_params['P2'],
                              var_params['P3'], var_params['P4'],
                              var_params['P5'], var_params['P6'],
                              var_params['P7'], var_params['P8'],
                              var_params['P9']])
        hkl = ['220', '311', '222',
               '004', '331', '224', '333',
               '511', '440']
        intensities = [var_params['I1'], var_params['I2'], var_params['I3'],
                       var_params['I4'], var_params['I5'], var_params['I6'],
                       var_params['I7'], var_params['I8'], var_params['I9']]
        return whole_pattern(self.q, var_params['D'],var_params['delta'],
                             8.385, positions, hkl, intensities,
                             var_params['eta'], var_params['sig'], self.bkg_poly)
                             #, var_params['bkgm'], var_params['bkgt'])

    def square_differences(self, params):
        calc = self.fit_function(params)
        res = sum((self.I - calc) ** 2) / sum(self.I ** 2)
        return np.array(res)

    def iteration_number(self, params, iter_n, resid):
        print(iter_n)

    def simulate(self):
        self.fit_curve = self.fit_function(self.params)

    def fit(self):
        result = minimize(self.square_differences, self.params, method='nelder-mead',
                          iter_cb=self.iteration_number)
        print(self.name)
        result.params.pretty_print()
        print("-----------------------------------------")
        self.result_params = result.params
        self.fit_curve = self.fit_function(self.result_params)


def A_S(L, D):
    return 1 - 3 / 2 * (L / D) + 1 / 2 * (L / D) ** 3


def A_APD(L, h, k, l, a0, delta, D):
    # unaffected peaks
    if (h + k) % 4 == 0 and (h + l) % 4 == 0 and (k + l) % 4 == 0:
        fac = 0

    # s.1 (220), (422)
    elif (h) % 2 == 0 and (k) % 2 == 0 and (l) % 2 == 0:
        fac = (2 * h + 2 * k + 2 * l) / (a0 * np.sqrt(2 * (h ** 2 + k ** 2 + l ** 2)))

    else:
        # s.2.1 (311)
        if (h + k) % 4 == 0 and (h + l) % 4 == 0:
            fac = (h + k + l) / (a0 * np.sqrt((h ** 2 + k ** 2 + l ** 2)))

        # s.2.2 (511), (333)
        else:
            fac = 4 * l / (a0 * np.sqrt(3 * (h ** 2 + k ** 2 + l ** 2)))

    return (1 - 2 * delta) ** (L * fac)


def A_IP(L, nu, sig):
    if sig == 0:
        return 1
    k = 1 / (1 + (0.677622) * (1 - nu) / nu)
    A = (1 - k) * np.exp(-np.pi ** 2 * sig ** 2 * L ** 2 / np.log(2))
    B = k * np.exp(-2 * np.pi * sig * L)
    return A + B


def profile_APB_size(Q, D, hkl, a0, delta, pos, nu, sig):
    h, k, l = hkl
    h, k, l = int(h), int(k), int(l)
    N = 1000
    L = np.linspace(0, D, N)
    L = L[:, np.newaxis]

    Fourier_coeff = A_S(L, D) * A_APD(L, h, k, l, a0, delta, D) * A_IP(L, nu, sig)

    exponent = ne.evaluate("exp(1j*L*(Q-pos)).real")
    I = ne.evaluate("Fourier_coeff*exponent")
    I = ne.evaluate("sum((I), axis=0)")
    return I / N * D


def whole_pattern(q, D, delta, a0, positions, hkl, intensities, nu, sig, bkg_poly):#, bkgm, bkgt):
    I = np.zeros_like(q)

    for hkl, pos, d, intensity in zip(hkl, positions, [delta] * len(hkl),
                                      intensities):
        I += intensity * profile_APB_size(q, D, hkl, a0, d, pos, nu, sig)

    return I + bkg_poly #(bkgm * q + bkgt)