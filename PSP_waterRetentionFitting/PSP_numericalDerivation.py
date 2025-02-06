# PSP_numericalDerivation.py
import numpy as np


def firstDerivative5Points(y):
    derivativeList = []
    for i in range(0, 2):
        derivativeList.append(np.nan)
    for i in range(2, len(y) - 2):
        dY = (1. / 12.) * (y[i - 2] - 8. * y[i - 1] + 8. * y[i + 1] - y[i + 2])
        derivativeList.append(dY)
    for i in range(len(y) - 2, len(y)):
        derivativeList.append(np.nan)
    return derivativeList
