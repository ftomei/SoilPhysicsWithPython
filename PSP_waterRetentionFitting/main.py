# PSP_waterRetentionFitting
from __future__ import print_function, division

import matplotlib.pyplot as plt
from PSP_readDataFile import readDataFile
from PSP_Marquardt import *
from PSP_numericalDerivation import *

def main():
    # read the experimental values
    myOutput, isFileOk = readDataFile("data/wheat/PAR.txt", 1, '\t', False)
    if not isFileOk:
        print('Wrong file: error reading row nr.', myOutput)
        return False
    waterPotential = myOutput[:, 0]
    waterContent = myOutput[:, 1]
    
    # select water retention curve
    print(CAMPBELL, ' Campbell')
    print(VAN_GENUCHTEN, ' van Genuchten')
    print(RESTRICTED_VG, ' van Genuchten with m = 1-1/n restriction')
    print(RESTRICTED_VG_BIMODAL, ' van Genuchten bimodal with m = 1-1/n restriction')
    print(IPPISCH_VG, ' Ippisch-van Genuchten')
    print(RESTRICTED_IPPISCH_VG, ' Ippisch-van Genuchten with m = 1-1/n restriction')
    print(CAMPBELL_IPPISCH_VG, ' Campbell-Ippisch-van Genuchten')

    waterRetentionCurve = 0
    while (waterRetentionCurve < CAMPBELL) or (waterRetentionCurve > LAST_MODEL_NR):
        waterRetentionCurve = float(input("Choose model type: "))
        if (waterRetentionCurve < CAMPBELL) or (waterRetentionCurve > LAST_MODEL_NR):
            print("wrong choice.")

    # initialize parameters
    thetaSatList = []
    previousWC = 0.0
    for i in range(len(waterContent)):
        wc = waterContent[i]
        if wc > previousWC:
            thetaSatList.append(wc)
        previousWC = wc

    minThetaS = min(thetaSatList)
    thetaS = sum(thetaSatList) / len(thetaSatList)
    thetaR = min(waterContent)

    # initial values
    air_entry = 1.0
    Campbell_b = 4.0
    VG_alpha = 1/air_entry
    VG_n = 1.2
    VG_m = 1. - 1./VG_n

    # bimodal initial parameters
    VG_alpha2 = VG_alpha * 0.9
    VG_n2 = VG_n * 1.1
    w = 0.5
      
    if waterRetentionCurve == CAMPBELL:
        b0 = np.array([thetaS, air_entry, Campbell_b], float)
        bmin = np.array([minThetaS, 0.1, 0.1], float)
        bmax = np.array([1., 20., 10.], float)
    elif waterRetentionCurve == VAN_GENUCHTEN:
        b0 = np.array([thetaS, thetaR, VG_alpha, VG_n, VG_m], float)
        bmin = np.array([minThetaS, 0., 0.001, 0.01, 0.01], float)
        bmax = np.array([1., thetaR, 10., 10., 1.], float)
    elif waterRetentionCurve == RESTRICTED_VG:
        b0 = np.array([thetaS, thetaR, VG_alpha, VG_n], float)
        bmin = np.array([minThetaS, 0., 0.001, 1.], float)
        bmax = np.array([1., thetaR, 10., 10.], float)
    elif waterRetentionCurve == RESTRICTED_VG_BIMODAL:
        b0 = np.array([thetaS, thetaR, VG_alpha, VG_alpha2, VG_n, VG_n2, w], float)
        bmin = np.array([minThetaS, 0., 0.001, 0.001, 1., 1., 0.], float)
        bmax = np.array([1., thetaR, 100., 100., 10., 10., 1.], float)
    elif waterRetentionCurve == IPPISCH_VG:
        b0 = np.array([thetaS, thetaR, air_entry, VG_alpha, VG_n, VG_m], float)
        bmin = np.array([minThetaS, 0., 0.1, 0.001, 0.01, 0.01], float)
        bmax = np.array([1., thetaR, 10., 10., 10., 1.], float)
    elif waterRetentionCurve == RESTRICTED_IPPISCH_VG:
        b0 = np.array([thetaS, thetaR, air_entry, VG_alpha, VG_n], float)
        bmin = np.array([minThetaS, 0., 0.1, 0.001, 1.], float)
        bmax = np.array([1., thetaR, 10., 10., 10.], float)
    elif waterRetentionCurve == CAMPBELL_IPPISCH_VG:
        b0 = np.array([thetaS, thetaR, air_entry, VG_alpha, VG_n], float)
        bmin = np.array([minThetaS, 0., 0.1, 0.01, 1.], float)
        bmax = np.array([1., thetaR, 10., 10., 10.], float)
    else:
        print("wrong choice.")
        return False

    # assign field capacity and wilting point
    isFieldCapacityCorrect = False
    fieldCapacity = np.zeros(1)     # [J kg-1]
    while not isFieldCapacityCorrect:
        fieldCapacity[0] = float(input("Insert field capacity value [J kg-1]: "))
        if (fieldCapacity[0] < 0) or (fieldCapacity[0] > 100):
            print('wrong choice.')
        else:
            isFieldCapacityCorrect = True

    wiltingPoint = np.zeros(1)
    wiltingPoint[0] = 1500          # [J kg-1]

    print("\nFitting")
    b = Marquardt(waterRetentionCurve, b0, bmin, bmax, waterPotential, waterContent)

    print("\nthetaS = ", b[0])
    if waterRetentionCurve == CAMPBELL:
        print("AirEntry = ", b[1])
        print("b = ", b[2])
    elif waterRetentionCurve == VAN_GENUCHTEN:
        print("thetaR = ", b[1])
        print("alpha = ", b[2])
        print("n = ", b[3])
        print("m = ", b[4])
    elif waterRetentionCurve == RESTRICTED_VG:
        print("thetaR = ", b[1])
        print("alpha = ", b[2])
        print("n = ", b[3])
    elif waterRetentionCurve == RESTRICTED_VG_BIMODAL:
        print("thetaR = ", b[1])
        print("alpha1 = ", b[2])
        print("alpha2 = ", b[3])
        print("n1 = ", b[4])
        print("n2 = ", b[5])
        print("w0 = ", b[6])
        print("w1 = ", 1 - b[6])
    elif waterRetentionCurve == IPPISCH_VG:
        print("thetaR = ", b[1])
        print("AirEntry = ", b[2])
        print("alpha = ", b[3])
        print("n = ", b[4])
        print("m = ", b[5])
    elif waterRetentionCurve == RESTRICTED_IPPISCH_VG:
        print("thetaR = ", b[1])
        print("AirEntry = ", b[2])
        print("alpha = ", b[3])
        print("n = ", b[4])
    elif waterRetentionCurve == CAMPBELL_IPPISCH_VG:
        print("thetaR = ", b[1])
        print("AirEntry = ", b[2])
        print("alpha = ", b[3])
        print("n = ", b[4])

    wcFieldCapacity = estimateRetention(waterRetentionCurve, b, fieldCapacity)
    wcWiltingPoint = estimateRetention(waterRetentionCurve, b, wiltingPoint)
    print()
    print("Field Capacity [m3 m-3] =", wcFieldCapacity)
    print("Wilting Point [m3 m-3] =", wcWiltingPoint)
    print("Plant Available Water [m3 m-3] =", wcFieldCapacity - wcWiltingPoint)

    plt.figure(figsize=(8, 6))

    # different colors for each series (maximum 4)
    colorList = ['r.', 'g.', 'y.', 'b.']
    colorIndex = 0
    previousWP = 0
    for i in range(len(waterPotential)):
        wp = waterPotential[i]
        if wp < previousWP and colorIndex < 3:
            colorIndex += 1
        plt.plot(waterPotential[i], waterContent[i], colorList[colorIndex])
        previousWP = wp

    plt.xscale('log')
    plt.xlabel('Water Potential [J kg$^{-1}$]', fontsize=14, labelpad=6)
    plt.xticks(size='16')
    plt.ylabel('Water Content [m$^{3}$ m$^{-3}$]', fontsize=14, labelpad=6)
    plt.tick_params(axis='both', which='major', labelsize=14, pad=6)
    plt.tick_params(axis='both', which='minor', labelsize=14, pad=6)

    # plot water retention curve
    myWP = np.logspace(-3, 8, 100)
    myWC = estimateRetention(waterRetentionCurve, b, myWP)
    plt.plot(myWP, myWC, 'k')

    plt.savefig('waterRetention.jpg')

    if waterRetentionCurve == RESTRICTED_VG_BIMODAL:
        # leggere file delle conducibilità
        k_sat = 0.01   # kg s m-3

        # plot conductivity
        myConductivity = estimateConductivity(waterRetentionCurve, b, k_sat, myWP)
        fig2, ax2 = plt.subplots()
        ax2.set_xlabel('Water Potential [J kg$^{-1}$]')
        ax2.set_ylabel('Hydraulic conductivity [kg s m$^{-3}$]')
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.plot(myWP, myConductivity, 'k')

         # plot soil pore radius
        degreeOfSaturation = getDegreeOfSaturation(waterRetentionCurve, b, np.flip(myWC))
        radius = getPoreRadius(np.flip(myWP), 30)
        degreeOfSaturationPdf = firstDerivative5Points(degreeOfSaturation)

        fig1, ax1 = plt.subplots()
        ax1.set_xlabel('radius (m)')
        ax1.set_xscale('log')
        ax1.set_ylabel('degree of saturation [-]')
        ax1.plot(radius, degreeOfSaturation, 'k')

        ax2 = ax1.twinx()           # instantiate a second axes that shares the same x-axis
        ax2.set_ylabel('pdf [-]')   # we already handled the x-label with ax1
        ax2.plot(radius, degreeOfSaturationPdf, 'r')

        fig1.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.savefig('pdf.jpg')

    plt.show()


main()    
