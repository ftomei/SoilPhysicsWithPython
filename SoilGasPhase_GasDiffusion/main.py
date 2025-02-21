#PSP_gasDiffusion
from math import exp
from PSP_ThomasAlgorithm import ThomasBoundaryCondition
import PSP_grid as grid
import matplotlib.pyplot as plt
import numpy as np

def gasSolver(boundaryLayerCond, boundaryOxygenConc, dg, 
                                respRate, totalDepth, n):
    a  = np.zeros(n+2, float)  
    b  = np.zeros(n+2, float)  
    c  = np.zeros(n+2, float)  
    d  = np.zeros(n+2, float) 
    g  = np.zeros(n+2, float) 
    u  = np.zeros(n+2, float)  
    co = np.zeros(n+2, float)  
    
    g[0] = boundaryLayerCond
    co[0] = boundaryOxygenConc
    # vector depth [m]
    z = grid.linear(n, totalDepth)
    
    # initialize matrix
    for i in range(1, n+1):
        u[i] = respRate * exp(-z[i] / 0.3) * (z[i + 1] - z[i - 1]) / 2.0
        if i < n:
            g[i] = dg / (z[i + 1] - z[i])
        else:
            g[i] = 0
        a[i + 1] = -g[i]
        b[i] = g[i - 1] + g[i]
        c[i] = -g[i]
        d[i] = u[i]

    d[1] = d[1] + g[0] * co[0]
    
    ThomasBoundaryCondition(a, b, c, d, co, 1, n)
    
    return(z, co)


def main():
    R = 8.3143                     
    n = 20                        
    totalDepth = 0.5               
    bulkDensity = 1200.            
    particleDensity = 2650.         
    waterContent = 0.2                 
    respRate = -0.001               
    oxygenDiff = 1.39e-5  #this is the value for carbon dioxide         
    temperature = 25.             
    atmPressure = 101.3           
    boundaryLayerCond = 0.01     
    
    # O2 concentration in air [g/m^3]
    boundaryOxygenConc = (0.21 * atmPressure * 1000. * 32. / 
                          (R * (temperature + 273.15))) 
    porosity = 1. - bulkDensity / particleDensity
    gasPorosity = porosity - waterContent
    
    # binary diffusion coefficient [m2/s]
    binaryDiffCoeff = (oxygenDiff * (101.3 / atmPressure) 
                * ((temperature + 273.15) / 273.15)**1.75)
    
    bg = 0.9           
    mg = 2.3           
    dg = binaryDiffCoeff * bg * gasPorosity**mg
    
    z, co = gasSolver(boundaryLayerCond, boundaryOxygenConc, 
                      dg, respRate, totalDepth, n)
      
    print ("node   depth [m]   Co [g\m^3]")
    for i in range(n + 2):
        print ("%3d    %6.2f      %.2f" %(i, z[i], co[i]))
    
    # plot results
    fig = plt.figure(figsize=(10,8))
    for i in range(n+1):
        plt.plot(co[i], -z[i], 'ko')
        
    plt.xlabel('Concentration [g m$^{-3}$]',fontsize=20,labelpad=8)
    plt.ylabel('Depth [m]',fontsize=20,labelpad=8)
    plt.tick_params(axis='both', which='major', labelsize=20,pad=8)
    plt.tick_params(axis='both', which='minor', labelsize=20,pad=8)
    plt.show()
main()
