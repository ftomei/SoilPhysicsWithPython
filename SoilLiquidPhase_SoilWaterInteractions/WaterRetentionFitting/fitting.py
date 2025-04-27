# PSP_waterRetentionFitting
from __future__ import print_function, division

import matplotlib.pyplot as plt
import numpy as np
import csv
import os
from PSP_readDataFile import readDataFile
from PSP_Marquardt import *
from PSP_numericalDerivation import *


def assignWeights(obsWaterPotential, obsWaterContent, userWeight):
    nrObsValues = len(userWeight)
    waterPotential = []
    waterContent = []

    # check user weights
    sumUserWeight = 0
    nrNoWeightData = 0
    for i in range(nrObsValues):
        if userWeight[i] > 0:
            sumUserWeight += userWeight[i]
        else:
            nrNoWeightData += 1
    if sumUserWeight > 1.0:
        print('Wrong weight sum: ', sumUserWeight)
        return False, waterPotential, waterContent

    # set weights
    weight = np.zeros(nrObsValues)
    minWeight = 1
    for i in range(nrObsValues):
        if userWeight[i] != 0:
            weight[i] = userWeight[i]
        else:
            weight[i] = (1.0 - sumUserWeight) / nrNoWeightData
        if weight[i] < minWeight:
            minWeight = weight[i]

    nrNewData = int(1.0 / minWeight)
    for i in range(nrObsValues):
        nrData = round(weight[i] * nrNewData)
        for n in range(nrData):
            waterPotential.append(obsWaterPotential[i])
            waterContent.append(obsWaterContent[i])

    return True, waterPotential, waterContent


def main():
    # reads the experimental values
    myOutput, isFileOk = readDataFile("data/bimodal_weight.txt", 1, '\t', False)
    if not isFileOk:
        nrWrongRow = myOutput + 1
        print('Wrong file: error reading row nr.', nrWrongRow)
        return False

    obsWaterPotential = myOutput[:, 0]
    obsWaterContent = myOutput[:, 1]

    # use the weights assigned by users
    if len(myOutput[0]) < 3:
        waterPotential = obsWaterPotential
        waterContent = obsWaterContent
    else:
        userWeight = myOutput[:, 2]
        isOk, waterPotential, waterContent = assignWeights(obsWaterPotential, obsWaterContent, userWeight)
        if not isOk:
            return False

    # select the water retention curve
    print(CAMPBELL, ' Campbell')
    print(VAN_GENUCHTEN, ' van Genuchten')
    print(RESTRICTED_VG, ' van Genuchten with m = 1-1/n restriction')
    print(RESTRICTED_VG_BIMODAL, ' van Genuchten bimodal with m = 1-1/n restriction')
    print(IPPISCH_VG, ' Ippisch-van Genuchten')
    print(RESTRICTED_IPPISCH_VG, ' Ippisch-van Genuchten with m = 1-1/n restriction')
    print(CAMPBELL_IPPISCH_VG, ' Campbell-Ippisch-van Genuchten')

    waterRetentionCurve = 0
    while (waterRetentionCurve < CAMPBELL) or (waterRetentionCurve > LAST_MODEL_NR):
        waterRetentionCurve = float(input("Choose water retention curve: "))
        if (waterRetentionCurve < CAMPBELL) or (waterRetentionCurve > LAST_MODEL_NR):
            print("wrong choice.")

    # initialize ThetaS and ThetaR from observed data
    thetaSatList = []
    previousWC = 0.0
    for i in range(len(obsWaterContent)):
        wc = obsWaterContent[i]
        if wc > previousWC:
            thetaSatList.append(wc)
        previousWC = wc

    minThetaS = min(thetaSatList)
    thetaS = sum(thetaSatList) / len(thetaSatList)
    thetaR = min(waterContent)

    # Initial values of the water retention curve parameters
    air_entry = 1.0
    Campbell_b = 4.0
    VG_alpha = 1 / air_entry
    VG_n = 1.2
    VG_m = 1. - 1./VG_n

    # Initial values of bimodal parameters
    VG_alpha2 = VG_alpha * 0.5
    VG_n2 = VG_n * 0.5
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
        bmin = np.array([minThetaS, 0., 0.0001, 0.0001, 1., 1., 0.], float)
        bmax = np.array([1., thetaR, 100., 100., 20., 20., 1.], float)
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

    # assigns the field capacity and the wilting point
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

    print("\nFitting...")
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

    # plot the observed water potential with different colors for each series (maximum 4)
    colorList = ['r.', 'g.', 'y.', 'b.']
    colorIndex = 0
    previousWP = 0
    for i in range(len(obsWaterPotential)):
        wp = obsWaterPotential[i]
        if wp < previousWP and colorIndex < 3:
            colorIndex += 1
        plt.plot(obsWaterPotential[i], obsWaterContent[i], colorList[colorIndex])
        previousWP = wp

    plt.xscale('log')
    plt.xlabel('Water Potential [J kg$^{-1}$]', fontsize=14, labelpad=6)
    plt.xticks(size='16')
    plt.ylabel('Water Content [m$^{3}$ m$^{-3}$]', fontsize=14, labelpad=6)
    plt.tick_params(axis='both', which='major', labelsize=14, pad=6)
    plt.tick_params(axis='both', which='minor', labelsize=14, pad=6)

    # plot the fitted water retention curve
    myWP = np.logspace(-3, 8, 200)
    myWC = estimateRetention(waterRetentionCurve, b, myWP)
    plt.plot(myWP, myWC, 'k')

    if not os.path.exists('output'):
        os.makedirs('output')
    plt.savefig('output/waterRetention.jpg')

    if waterRetentionCurve == RESTRICTED_VG_BIMODAL:
        # read water conductivity
        myInput, isFileOk = readDataFile("data/bimodal_K.csv", 1, ',', False)
        if not isFileOk:
            print('Wrong file: error reading row nr.', myInput)
            return False
        waterPotential_k = myInput[:, 0]           # [J kg-1]
        waterConductivity = myInput[:, 1]          # [kg s m-3]
        k_sat = max(waterConductivity)

        # plot the estimated water conductivity
        myConductivity = estimateConductivity(waterRetentionCurve, b, k_sat, myWP)
        f2, fig2 = plt.subplots()
        fig2.set_xlabel('Water Potential [J kg$^{-1}$]')
        fig2.set_ylabel('Hydraulic conductivity [kg s m$^{-3}$]')
        fig2.set_xscale('log')
        fig2.set_yscale('log')
        fig2.plot(myWP, myConductivity, 'k')

        # save the estimated data to csv file
        outputFilename = "output/fitting.csv"
        with open(outputFilename, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["water potential", "water content", "water conductivity"])
            for c1, c2, c3 in zip(myWP, myWC, myConductivity):
                writer.writerow([c1, c2, c3])
        print('\nOutput file:', outputFilename, '\n')

        # plot the observed conductivity with different colors for each series (maximum 4)
        colorList = ['r.', 'g.', 'y.', 'b.']
        colorIndex = 0
        previousWP = 0
        for i in range(len(waterPotential_k)):
            wp = waterPotential_k[i]
            if wp < previousWP and colorIndex < 3:
                colorIndex += 1
            fig2.plot(waterPotential_k[i], waterConductivity[i], colorList[colorIndex])
            previousWP = wp
        plt.savefig('output/waterConductivity.jpg')

        # plot soil pore radius and pdf
        degreeOfSaturation = getDegreeOfSaturation(waterRetentionCurve, b, np.flip(myWC))
        radius = getPoreRadius(np.flip(myWP), 30)
        degreeOfSaturationPdf = firstDerivative5Points(degreeOfSaturation)

        # search for relative maximum values of pdf
        secondDerivative = firstDerivative5Points(degreeOfSaturationPdf)
        isFirst = True
        previousSign = 1
        for i in range(0, len(secondDerivative)):
            if secondDerivative != math.nan:
                if secondDerivative[i] < 0:
                    currentSign = -1
                else:
                    currentSign = 1

                if not isFirst:
                    if currentSign == -1 and previousSign == 1:
                        print("Relative maximum of Pdf: ", radius[i])
                else:
                    isFirst = False
                previousSign = currentSign

        f3, fig3 = plt.subplots()
        fig3.set_xlabel('radius (m)')
        fig3.set_xscale('log')
        fig3.set_ylabel('degree of saturation [-]')
        fig3.plot(radius, degreeOfSaturation, 'k')

        fig4 = fig3.twinx()             # instantiate a second axes that shares the same x-axis
        fig4.set_ylabel('pdf [-]')      # we already handled the x-label with ax1
        fig4.plot(radius, degreeOfSaturationPdf, 'r')

        f3.tight_layout()               # otherwise, the right y-label is slightly clipped
        plt.savefig('output/pdf.jpg')

    plt.show()


main()    
