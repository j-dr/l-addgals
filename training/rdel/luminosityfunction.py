from __future__ import print_function, division
from copy import copy
import numpy as np


class LuminosityFunction(object):

    def __init__(self, params):

        self.params = params
        self.lf = None

    def genLuminosityFunction(self, lums, zs):

        self.lf = np.zeros((len(lums), len(zs)))

        for i, z in enumerate(zs):
            zp = self.evolveParams(z)
            self.lf[:,i] = self.calcNumberDensity(zp, lums)

    def genLuminosityFunctionZ(self, lums, z):
        zp = self.evolveParams(z)
        return self.calcNumberDensity(zp, lums)

    def calcNumberDensity(self, p, lums):
        pass

    def evolveParams(self, z):
        pass


class DSGLuminosityFunction(LuminosityFunction):

    def __init__(self, params=None):

        if params is None:
            params = np.zeros(8)
            mstar = -20.44
            mstar0 = -20.310

            params[0] = 0.0156  #phistar1
            params[1] = -0.166  #alpha1
            params[2] = 0.00671 #phistar2
            params[3] = -1.523  #alpha2
            params[4] = -19.88  #mstar
            params[5] = 3.08e-5 #phistar3
            params[6] = -21.72  #M_hi
            params[7] = 0.484   #sigma_hi

            mr_shift = mstar - mstar0
            params[4] += mr_shift
            params[6] += mr_shift

        LuminosityFunction.__init__(self,params)

    def evolveParams(self, z):
        zp = copy(self.params)

        phistar = 10 ** (-1.79574 + (-0.266409 * z))
        phistar_rat = phistar/self.params[0]

        zp[0] *= phistar_rat
        zp[2] *= phistar_rat
        zp[5] *= phistar_rat

        return zp

    def calcNumberDensity(p, lums):
        """
        Sum of a double schechter function and a gaussian.
        m -- magnitudes at which to calculate the number density
        p -- Function parameters. Order
             should be phi^{star}_{1}, M^{star}, \alpha_{1},
             phi^{star}_{2}, M^{star}, \alpha_{2}, \phi_{gauss},
             \M_{gauss}, \sigma_{gauss}
        """
        phi = 0.4 * np.log(10) * np.exp(-10**(-0.4 * (lums - p[4]))) * \
            (p[0] * 10 ** (-0.4 * (lums - p[4])*(p[1]+1)) + \
            p[2] * 10 ** (-0.4 * (lums - p[4])*(p[3]+1))) + \
            p[5] / np.sqrt(2 * np.pi * p[7] ** 2) * \
            np.exp(-(lums - p[6]) ** 2 / (2 * p[7] ** 2))

        return phi
