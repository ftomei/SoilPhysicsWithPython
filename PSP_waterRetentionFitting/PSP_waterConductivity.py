# PSP_waterConductivity.py

from __future__ import division
import math
import numpy as np
from PSP_waterRetention import *


def VanGenuchtenRestrictedBimodalConductivity(v, k_sat, psi, cond):
    alpha1 = v[2]
    alpha2 = v[3]

    n1 = v[4]
    n2 = v[5]
    m1 = 1. - (1. / n1)
    m2 = 1. - (1. / n2)

    w1 = v[6]
    w2 = 1 - w1

    for i in range(len(psi)):
        if psi[i] <= 0:
            k = k_sat
        else:
            Se1 = (1. + (alpha1 * abs(psi[i])) ** n1) ** (-m1)
            Se2 = (1. + (alpha2 * abs(psi[i])) ** n2) ** (-m2)
            num1 = alpha1 * (1 - (1 - Se1 ** (1 / m1)) ** m1)
            num2 = alpha2 * (1 - (1 - Se2 ** (1 / m2)) ** m2)
            num = w1 * num1 + w2 * num2
            den = w1 * alpha1 + w2 * alpha2
            Se = w1 * Se1 + w2 * Se2
            k = k_sat * Se ** 0.5 * (num / den) ** 2
        cond[i] = k


def VanGenuchtenRestrictedConductivity(v, k_sat, psi, cond):
    alpha = v[2]
    n = v[3]
    m = 1. - (1. / n)

    for i in range(len(psi)):
        if psi[i] <= 0:
            k = k_sat
        else:
            Se = (1. + (alpha * abs(psi[i])) ** n) ** (-m)
            num = alpha * (1 - (1 - Se ** (1 / m)) ** m)
            k = k_sat * Se ** 0.5 * (num / alpha) ** 2
        cond[i] = k
