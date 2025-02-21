# PSP_basicProperties.py

def computePorosity(bulkDensity, massWetness):
    waterDensity = 1000.
    particleDensity = 2650.
    porosity = 1 - (bulkDensity / particleDensity)
    voidRatio = porosity / (1 - porosity)
    waterContent = massWetness * (bulkDensity / waterDensity)
    gasPorosity = porosity - waterContent
    degreeSaturation = waterContent / porosity
    print("\nTotal porosity [m^3/m^3] = ", format(porosity, '.3f'))
    print("Void ratio [m^3/m^3] = ", format(voidRatio, '.3f'))
    print("Volumetric water content [m^3/m^3] = ", format(waterContent, '.3f'))
    print("Gas filled porosity [m^3/m^3] = ", format(gasPorosity, '.3f'))
    print("Degree of saturation [-] =", format(degreeSaturation, '.3f'))
    return


def computeSaturationWetness(bulkDensity):
    waterDensity = 1000
    particleDensity = 2650
    porosity = 1 - (bulkDensity / particleDensity)
    return porosity / (bulkDensity / waterDensity)


def main():
    bulkDensity = float(input("Bulk Density [kg/m^3] = "))
    massWetness = float(input("Mass wetness [kg/kg] = "))

    satMassWetness = computeSaturationWetness(bulkDensity)
    if (massWetness >= 0) and (massWetness < satMassWetness):
        computePorosity(bulkDensity, massWetness)
    else:
        print("Wrong mass wetness! value at saturation = ", satMassWetness)

    input("Press a key to continue..")


main()
