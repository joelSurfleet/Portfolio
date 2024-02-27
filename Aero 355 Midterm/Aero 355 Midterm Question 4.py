# Aero 355 Midterm Question 4
# Joel Surfleet

import matplotlib.pyplot as plt

k = 1.38e-23 # J/K
g = 9.81 # m/s2
R = 8.314 # J/(mol*K)
pi = 3.14159265358979

f = open(r"C:\Users\joels\Documents\Python\Aero 355 Midterm\AtmosData.txt")
data = f.read()
f.close()
data = data.splitlines()

h = [None] * len(data)  # Height [km]

molMassMin = [None] * len(data) # molar mass @ solar min [kg/mol]
molMassMax = [None] * len(data) # molar mass @ solar max [kg/mol]

crossSectMin = [None] * len(data) # cross sectional area @ solar min [m^2]
crossSectMax = [None] * len(data) # cross sectional area @ solar max [m^2]

for i in range(len(data)):
    h[i] = float(data[i].split()[0])

    molMassMin[i]   = float(data[i].split()[1])
    molMassMax[i]   = float(data[i].split()[2])
    
    crossSectMin[i] = float(data[i].split()[3])
    crossSectMax[i] = float(data[i].split()[4])

T0 = 180 # K
P0 = 1 #Pa
h0 = h[0]

lapseMin = (700 - 180) / (400 - 80)
lapseMax = (1800 - 180) / (400 - 80)

tempMin = [None] * len(data)
tempMax = [None] * len(data)

pressMin = [None] * len(data)
pressMax = [None] * len(data)

densityMin = [None] * len(data)
densityMax = [None] * len(data)

thermVelocMin = [None] * len(data)
thermVelocMax = [None] * len(data)

mfpMin = [None] * len(data)
mfpMax = [None] * len(data)

def lapsePressure(L,h,M):
    return P0 * (1 + L/T0 * (h - h0)) ** (-(M*g)/(R*L))

def numberDensity(P,T):
    return P/(k*T)

def thermalVelocity(T,M):
    return ((8*k*T)/(pi*M)) ** (1/2)

def meanFreePath(density,r):
    sigma = pi*r**2
    return 1/(4*density*sigma)

for i in range(len(data)):
    tempMin[i] = T0 + lapseMin * (h[i] - h[0]) # K
    tempMax[i] = T0 + lapseMax * (h[i] - h[0]) # K

    pressMin[i] = lapsePressure(lapseMin,h[i],molMassMin[i]) # Pa
    pressMax[i] = lapsePressure(lapseMax,h[i],molMassMax[i]) # Pa

    densityMin[i] = numberDensity(pressMin[i],tempMin[i])
    densityMax[i] = numberDensity(pressMax[i],tempMax[i])

    thermVelocMin[i] = thermalVelocity(tempMin[i],molMassMin[i])
    thermVelocMax[i] = thermalVelocity(tempMax[i],molMassMax[i])
    
    mfpMin[i] = meanFreePath(densityMin[i],crossSectMin[i])
    mfpMax[i] = meanFreePath(densityMax[i],crossSectMax[i])

plt.figure()
plt.plot(h,pressMin,label="Solar Min")
plt.plot(h,pressMax,label="Solar Max")
plt.xlabel("Height [km]")
plt.ylabel("Pressure [Pa]")
plt.title("Pressure vs Altitude\n@ 400 km, Solar Min: " + str(pressMin[-1]) +" Pa\n@ 400 km, Solar Max: " + str(pressMax[-1]) + " Pa")
plt.legend()
plt.grid()

plt.figure(2)
plt.plot(h,densityMin,label="Solar Min")
plt.plot(h,densityMax,label="Solar Max")
plt.xlabel("Height [km]")
plt.ylabel("Density [1/m^3]")
plt.title("Density vs Altitude\n@ 400 km, Solar Min: " + str(densityMin[-1]) +" 1/m^3\nSolar Max: " + str(densityMax[-1]) + " 1/m^3")
plt.legend()
plt.grid()

plt.figure(3)
plt.plot(h,thermVelocMin,label="Solar Min")
plt.plot(h,thermVelocMax,label="Solar Max")
plt.xlabel("Height [km]")
plt.ylabel("Thermal Velocity [m/s]")
plt.title("Thermal Velocity vs Altitude\n@ 400 km, Solar Min: " + str(thermVelocMin[-1]) +" m/s\nSolar Max: " + str(thermVelocMax[-1]) + " m/s")
plt.legend()
plt.grid()

plt.figure(4)
plt.plot(h,mfpMin,label="Solar Min")
plt.plot(h,mfpMax,label="Solar Max")
plt.xlabel("Height [km]")
plt.ylabel("Mean Free Path [m]")
plt.title("Mean Free Path vs Altitude\n@ 400 km, Solar Min: " + str(mfpMin[-1]) +" m\nSolar Max: " + str(mfpMax[-1]) + " m")
plt.legend()
plt.grid()

plt.show()