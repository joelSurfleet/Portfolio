# Aero 355 midterm Question 7
# Joel Surfleet

import numpy as np
import matplotlib.pyplot as plt

# m = 10 # kg

TML = 0.25/100 # %
CVCM = 0.01/100 # %

h = 6.62607015e-34 # J.s
k = 1.38e-23 # J/K

tau0 = 3e-12 # s

# a)
def residenceTime(activationEnergy,T):
    return tau0*np.exp(activationEnergy/(R*T)) # s

T0 = 50
T = [T0 + 10*x for x in range(int((300-50)/10 + 1))]

R = 1.9858775e-3

E5  = [None] * len(T)
E10 = [None] * len(T)
E15 = [None] * len(T)
E20 = [None] * len(T)

for i in range(len(T)):
    E5[i]  = residenceTime(5, T[i])
    E10[i] = residenceTime(10,T[i])
    E15[i] = residenceTime(15,T[i])
    E20[i] = residenceTime(20,T[i])

plt.figure(1)
plt.semilogy(T,E5, label="5 kcal/mole")
plt.semilogy(T,E10,label="10 kcal/mole")
plt.semilogy(T,E15,label="15 kcal/mole")
plt.semilogy(T,E20,label="20 kcal/mole")
plt.xlabel("Temperature [K]")
plt.ylabel("Residence Time [s]")
plt.title("Residence Time vs Temperature\nResidence Time @ 300K, 20 kcal/mol: " + str(E20[-1]) + " s")
plt.legend()
plt.grid()

A = 0.1 # m^2
m = 5 # g

E = 15

q0 = (TML*m)/(2*np.exp(-15/(R*298))) # g/day

def massLoss(t2,t1,T,activationEnergy):
    return 2*q0*np.exp(-activationEnergy/(R*T))*(t2**(1/2) - t1**(1/2)) # g/day

rho = 1 # g/cm^3

def arrivalRate(viewFactor,mDot):
    return viewFactor * mDot/rho

num = 20
time = range(8)
massLoss0 = [None] * len(time)
massLoss1 = [None] * len(time)
massLoss2 = [None] * len(time)
arrRate = [[0 for i in range(8)] for j in range(num)]

slope = (0.03-0.0003) / (num-1)

for i in range(1,len(time)):
    massLoss0[i] = massLoss(time[i],time[i-1],273,15)
    massLoss1[i] = massLoss(time[i],time[i-1],323,15)
    massLoss2[i] = massLoss(time[i],time[i-1],373,15)

    mdot250k = massLoss(time[i],time[i-1],250,15)

    for j in range(num):
        arrRate[j][i] = arrivalRate(0.0003+slope*(j),mdot250k)
        arrRate[j][0] = None

plt.figure(2)
plt.semilogy(time,massLoss0,label="0 degrees C")
plt.semilogy(time,massLoss1,label="50 degrees C")
plt.semilogy(time,massLoss2,label="100 degrees C")
plt.xlabel("Time [days]")
plt.ylabel("Daily Mass Loss [g/day]")
plt.title("Daily Mass Loss vs Time\nMass Loss @ 100C, 7th day: " + str(massLoss2[-1]) + " g")
plt.legend()
plt.grid()

plt.figure(3)
for i in range(num):
    plt.semilogy(time,arrRate[i],label=str(0.0003+slope*(i)))
plt.xlabel("Time [days]")
plt.ylabel("Daily Arrival Rate [cm^3/day]")
plt.title("Daily Arrival rate vs Time for various View Factors\nArrival Rate @ 100C, 7th day: " + str(arrRate[-1][-1]) + " cm^3/day")
plt.legend(loc=1)
plt.grid()

plt.show()