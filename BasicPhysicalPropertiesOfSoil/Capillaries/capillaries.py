from vpython import *
from math import sqrt
from numpy import *
scene = canvas(width = 1000, height = 800, exit=True)
scene.background = color.white

n = 250
x = zeros(n)
z = zeros(n)
r = zeros(n)
for i in range(n):
 check = 0
 while(check == 0):
  x[i] = random.random_sample()-0.5
  z[i] = random.random_sample()-0.5
  check = 1
  for j in range(i):
    if ((x[i]-x[j])**2.+(z[i]-z[j])**2. < r[j]**2.):
      check = 0
  if not(x[i]**2.+z[i]**2.<0.2**2): check = 0
  r[i] = ((random.random_sample()*4-2)**2.+0.5)/4*0.1
  for j in range(i):
    if ((x[i]-x[j])**2.+(z[i]-z[j])**2. < (r[j]+r[i])**2.):
       r[i] = sqrt((x[i]-x[j])**2.+(z[i]-z[j])**2.)-r[j]
for i in range(n):
  cyl=cylinder(pos=vector(x[i],-1.5,z[i]), axis=vector(0,2,0), radius=r[i])
