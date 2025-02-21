#PSP_sedimentation
from __future__ import print_function, division
from vpython import *
from PSP_readDataFile import readDataFile

particleDensity = 2580
g = 9.80665
micro = 0.000001

def getSedimentationTime(solutionDensity, solutionViscosity, particleDiameter, z):
    numer = 18 * solutionViscosity * z
    denom = (g * (particleDensity - solutionDensity) * (particleDiameter * micro)**2)
    return (numer / denom)

# [m]
def getSedimentationDepth(solutionDensity, solutionViscosity, particleDiameter, time):
    numer = (time * g * (particleDensity - solutionDensity) * (particleDiameter * micro)**2)
    denom = 18 * solutionViscosity
    return (numer / denom)  

# [kg m^-3] 
def getSolutionDensity(liquidDensity, concentration):                       
	return (liquidDensity * (1 + 0.63*(concentration / 1000)))

# [kg s^-1 m^-1]
def getSolutionViscosity(liquidViscosity, concentration):                   
	return (liquidViscosity * (1 + 4.25*(concentration / 1000)))                               

def initializeScene(heightCylinder, radiusParticle):
    scene = canvas(width = 1000, height = 800, exit=True)
    scene.background = color.white
    scene.center = vector(0, 0, -heightCylinder/2)
    scene.forward = vector(0, 1, -0.5)
    
    cyl = cylinder(pos = vector(0,0,0), radius = 0.05)
    cyl.axis = vector(0,0,-heightCylinder)                                    
    cyl.color = color.red
    cyl.opacity = 0.2

    particle = sphere(pos = vector(0,0,0), color = color.black)
    particle.radius = min(radiusParticle * 0.001, 0.02) 
    return(particle)
    
    
def main():
    A, isFileOk = readDataFile("waterDensity.dat", 2, '\t', False)
    if ((not isFileOk) or (len(A[0]) < 3)):
        print("Wrong file!\nMissing data or wrong delimiter in line:", A+1)
        return()
        
    temperature = A[:,0]  
    density = A[:,1]
    viscosity = A[:,2]
    print ("Available temperature:", temperature)
    
    isGoodChoice = False
    while not isGoodChoice:
        t = float(input("\nInput temperature: "))
        for i in range(len(temperature)):
            if (temperature[i] == t):
                isGoodChoice = True
                liquidDensity = density[i]
                liquidViscosity = viscosity[i]
        if not isGoodChoice: 
            print ("Warning: not available value")
    
    concentration = float(input("Input concentration [g/l]: "))
    diameter = float(input("Input particle diameter [micrometers]: "))

    solutionDensity = getSolutionDensity(liquidDensity, concentration)
    solutionViscosity = getSolutionViscosity(liquidViscosity, concentration)
    
# [m]
    heightCylinder = 0.3    
    particle = initializeScene(heightCylinder, diameter/2.0)
            
# [s]
    time = 0                
# [m]
    depth = 0               
    while (depth < heightCylinder):
# wait 1/100 of second
        rate(100)    
        depth = getSedimentationDepth(solutionDensity, solutionViscosity, diameter, time)
        particle.pos.z = -depth
        time += 1
        
    print ("seconds:", time)    
main()
