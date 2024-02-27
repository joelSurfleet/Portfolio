# Aero 355 midterm Question 6
# Joel Surfleet

import matplotlib.pyplot as plt

f = open(r"C:\Users\joels\Documents\Python\Aero 355 Midterm\AODensity.txt")
data = f.read()
f.close()
data = data.splitlines()

for i in range(len(data)):
    h = float(data[i].split()[0])
    if h == 320:
        ngMin = float(data[i].split()[1])
        ngMax = float(data[i].split()[2])
        break

depth = 0.04 # cm
E = 10.5e-24  # Erosion yield in cm^3
E = E / (100 ** 3) # Erosion yield in m^3

def circularOrbitSpeed(z):
    a = z + 6378 # km
    return (398600 / a) ** (1/2) # km/s

V = 1000*circularOrbitSpeed(h) # m/s

def erosionRate(ng):
    dxdt = E * ng * V # m/s
    dxdt *= 86400 * 100 # cm/day
    return dxdt

time = range(366)

erosionMin = [erosionRate(ngMin) * day for day in time] # cm/day
erosionMax = [erosionRate(ngMax) * day for day in time] # cm/day

plt.figure(1)
plt.plot(time,erosionMin,label="Solar Min")
plt.plot(time,erosionMax,label="Solar Max")
plt.xlabel("Time [days]")
plt.ylabel("Erosion [cm]")
plt.title("Erosion vs Time\nErosion after 1 year @ Solar Min: " + str(erosionMin[-1]) + " cm\nErosion after 1 year @ Solar Max: " + str(erosionMax[-1]) + " cm")
plt.legend()
plt.grid()
plt.show()

timeMin = 0.04 / erosionRate(ngMin) # days
timeMax = 0.04 / erosionRate(ngMax) # days

print("Days until interconnect is eroded @ Solar Min:",timeMin)
print("Days until interconnect is eroded @ Solar Max:",timeMax)