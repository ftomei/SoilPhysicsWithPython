# PSP_waterRetention.py
from __future__ import division
import math
import numpy as np

CAMPBELL = 1
VAN_GENUCHTEN = 2
RESTRICTED_VG = 3
RESTRICTED_VG_BIMODAL = 4
IPPISCH_VG = 5
RESTRICTED_IPPISCH_VG = 6
CAMPBELL_IPPISCH_VG = 7
LAST_MODEL_NR = CAMPBELL_IPPISCH_VG

DEG_TO_RAD = 0.0174532925


def Campbell(v, psi, theta):
    thetaS = v[0]
    he = v[1]
    Campbell_b = v[2]
    for i in range(len(psi)):
        if psi[i] <= he:
            theta[i] = thetaS
        else:
            Se = (psi[i] / he) ** (-1. / Campbell_b)
            theta[i] = Se * thetaS


def VanGenuchten(v, psi, theta):
    thetaS = v[0]
    VG_thetaR = v[1]
    VG_alpha = v[2]
    VG_n = v[3]
    VG_m = v[4]
    for i in range(len(psi)):
        Se = 1. / pow(1. + pow(VG_alpha * psi[i], VG_n), VG_m)
        theta[i] = Se * (thetaS - VG_thetaR) + VG_thetaR


def VanGenuchtenRestricted(v, psi, theta):
    thetaS = v[0]
    VG_thetaR = v[1]
    VG_alpha = v[2]
    VG_n = v[3]
    VG_m = 1. - (1. / VG_n)
    for i in range(len(psi)):
        if psi[i] <= 0:
            Se = 0
        else:
            Se = (1. + (VG_alpha * abs(psi[i])) ** VG_n) ** (-VG_m)
        theta[i] = Se * (thetaS - VG_thetaR) + VG_thetaR


def VanGenuchtenRestrictedBimodal(v, psi, theta):
    thetaS = v[0]
    VG_thetaR = v[1]
    VG_alpha1 = v[2]
    VG_alpha2 = v[3]

    VG_n1 = v[4]
    VG_n2 = v[5]
    VG_m1 = 1. - (1. / VG_n1)
    VG_m2 = 1. - (1. / VG_n2)

    # normalized weight
    # v[6] = v[6] / (v[6] + v[7])
    # v[7] = 1 - v[6]
    # w2 = v[7]
    w1 = v[6]
    w2 = 1 - w1

    for i in range(len(psi)):
        if psi[i] <= 0:
            Se = 0
        else:
            Se1 = (1. + (VG_alpha1 * abs(psi[i])) ** VG_n1) ** (-VG_m1)
            Se2 = (1. + (VG_alpha2 * abs(psi[i])) ** VG_n2) ** (-VG_m2)
            Se = w1 * Se1 + w2 * Se2
        theta[i] = Se * (thetaS - VG_thetaR) + VG_thetaR


def IppischVanGenuchten(v, psi, theta):
    thetaS = v[0]
    VG_thetaR = v[1]
    he = v[2]
    VG_alpha = v[3]
    VG_n = v[4]
    VG_m = v[5]
    VG_Sc = (1. + (VG_alpha * he) ** VG_n) ** VG_m
    for i in range(len(psi)):
        if psi[i] <= he:
            Se = 1.0
        else:
            Se = VG_Sc * (1. + (VG_alpha * abs(psi[i])) ** VG_n) ** (-VG_m)
        theta[i] = Se * (thetaS - VG_thetaR) + VG_thetaR


def RestrictedIppischVanGenuchten(v, psi, theta):
    thetaS = v[0]
    VG_thetaR = v[1]
    he = v[2]
    VG_alpha = v[3]
    VG_n = v[4]
    VG_m = 1. - (1. / VG_n)
    VG_Sc = (1. + (VG_alpha * he) ** VG_n) ** VG_m
    for i in range(len(psi)):
        if psi[i] <= he:
            Se = 1.0
        else:
            Se = VG_Sc * (1. + (VG_alpha * abs(psi[i])) ** VG_n) ** (-VG_m)
        theta[i] = Se * (thetaS - VG_thetaR) + VG_thetaR


def CampbellIppischVanGenuchten(v, psi, theta):
    thetaS = v[0]
    VG_thetaR = v[1]
    he = v[2]
    VG_alpha = v[3]
    VG_n = v[4]
    VG_m = 1. - (1. / VG_n)
    VG_Sc = (1. + (VG_alpha * he) ** VG_n) ** VG_m
    for i in range(len(psi)):
        if psi[i] <= he:
            Se = 1.0
        else:
            Se = VG_Sc * (1. + (VG_alpha * abs(psi[i])) ** VG_n) ** (-VG_m)
        residual = VG_thetaR * (1 - (np.log(VG_alpha * psi[i] + 1.0) / np.log(VG_alpha * (10 ** 6) + 1.0)))
        theta[i] = max(0.0, Se * (thetaS - residual) + residual)


def getDegreeOfSaturation(waterRetentionCurve, parameters, waterContent):
    degOfSat = np.zeros(len(waterContent))

    if waterRetentionCurve == CAMPBELL:
        thetaS = parameters[0]
        for i in range(len(waterContent)):
            degOfSat[i] = waterContent[i] / thetaS
    else:
        thetaS = parameters[0]
        thetaR = parameters[1]
        for i in range(len(waterContent)):
            degOfSat[i] = (waterContent[i] - thetaR) / (thetaS - thetaR)

    return degOfSat


# [m]
def getPoreRadius(psi, angle):
    radius = np.zeros(len(psi))
    waterSurfaceTension = 72. * 0.001       # [N m-1]
    cosBeta = math.cos(angle * DEG_TO_RAD)
    waterDensity = 1000.                    # [kg m-3]

    for i in range(len(psi)):
        radius[i] = (2. * waterSurfaceTension * cosBeta) / (waterDensity * psi[i])

    return radius
