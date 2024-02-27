# Aero 355 midterm Question 5
# Joel Surfleet

import matplotlib.pyplot as plt

k = 1.38e-23 # J/K
g = 9.81 # m/s^2
R = 8.314 # J/(ml*K)

# a)
P0 = 1 # Pa
T0 = 273 - 80 # K
L = 3 # K/km

molMass = 16 # g/mol

h0 = 80 # km

h = [h0 + 20 * x for x in range(37)]

T = [T0 + L * (x - h0) for x in h]

rho0 = P0 / (R * T0)

def lapseDensity(h):
    return rho0 * (1 + L/T0 * (h - h0)) ** (-(molMass*g)/(R*L))

rho = [lapseDensity(x) for x in h]

plt.figure(1)
plt.semilogy(h,rho)
plt.xlabel("Height [km]")
plt.ylabel("Density [kg/m^3]")
plt.title("Density vs Altitude\nDensity @ 800 km: " + str(rho[-1]))
plt.grid()

# b)

areaSpacecraft = 100 # cm^2
areaSpacecraft = 100 / (100**2) # m^2
m = 1 # kg
Cd = 2.2 # Assumed to be 2.2 if not given

mu = 398600
rEarth = 6378

size = int(((10 - 0.01) / 0.01) + 1)

areaTotal = [areaSpacecraft + 0.01 * x for x in range(size)]

value = -rho[-1] * (Cd / m) * (mu * (rEarth + h[-1])) ** (1/2)

def dadt(A):
    return -rho[-1] * (Cd * A / m) * (mu * (rEarth + h[-1])) ** (1/2)

H = (R * T[-1])/(molMass * g)

time = [-H / dadt(A) for A in areaTotal]

time = [x / 86400 for x in time]

print(areaTotal[0])

plt.figure(2)
plt.semilogy(areaTotal,time)
plt.xlabel("Area [M^2]")
plt.ylabel("Time [days]")
plt.title("Lifetime vs Area\nLifetime of spacecraft (A = 100 cm^2): " + str(time[0]) + " days\nLifetime with full solar sail (A = "+ str(int(areaTotal[-1])) + " m^2): " + str(time[-1]) + " days")
plt.grid()

plt.show()