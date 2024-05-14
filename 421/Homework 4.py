from Homework4Functions import *
import matplotlib.pyplot as plt

# TODO
# - Fis Part 3 X-Scale
# - Do RMS on results from part 4, use that to find 90% confidence interval
# - Find largest square you can fit in the 90% confidence interval

err = {
    'time': 0.01,   # s     Local Sidereal Time
    'lati': 1e-4,   # lat   Target Location Latitude
    'long': 1e-4,   # lon   Target Location Longitude
    'elev': 10,     # mz    Target Location Elevation
    'loca': 3,      # mxyz  Spacecraft Location Knowledge (in/cross-track,radial)
    'V_xy': 0.002,  # m/sxy Spacecraft Velocity Knowledge (in/cross-track)
    'V__z': 0.007,  # m/sz  Spacecraft Velocity Knowledge (radial)
    'sPos': 0.01,   # mxy   Sensor Mounting Location
    'sAng': 1e-4,   # deg   Sensor Mounting Angle
}

# convert m to km
err['elev'] = err['elev'] / 1000
err['loca'] = err['loca'] / 1000
err['V_xy'] = err['V_xy'] / 1000
err['V__z'] = err['V__z'] / 1000
err['sPos'] = err['sPos'] / 1000

JD = 2458981.1666667 # Julian Date

r0 = [5980.8297,-1282.3184,4125.8019] # km
v0 = [1.540887,7.186813,0] # km/s

T0 = [35.3,-120.8,0.2] # Lat Long Alt

n = 1000 # number of iterations

g,gMag,_ = geoMonteSim(r0, v0, T0, JD, err, n) # run the sim!

# Answer to Part 1
print("Part 1")

# Take standard deviation, use that to find the 90% confidence interval in meters
gSTD = std(gMag) # km
g90 = 1.645*gSTD*1000 # m

print(g90)

gRad = g90 # store this for part 6

# Answer to Part 2
print("Part 2")

mu = 398600 # need this to recalc velocity

# Find rhat so we can vary its magnitude
rMag = norm(r0)
rHat = [r0[i] / rMag for i in range(3)]

# initialize stuff we need to store
g90 = [0 for i in range(32)]
pMag = [0 for i in range(32)]
rMag = [0 for i in range(32)]

# R_EARTH
R_EARTH = 6378

for i in range(32):
    # store values of r to step through so you can graph them later
    rMag[i] = R_EARTH + 1000 + i * 250

    # use magnitude to get an actual vector
    r = [rHat[j] * rMag[i] for j in range(3)]

    # I don't think we actually need this guy but whatever
    v = [(mu/r[j])**(1/2) for j in range(3)]

    # SIMULATE
    g,gMag,Tg = geoMonteSim(r, v, T0, JD, err, n)

    # Recalc p vector for graphing
    p = [r[i] - Tg[i] for i in range(3)]
    pMag[i] = norm(p)

    # Take standard deviation, use that to find the 90% confidence interval in meters
    gSTD = std(gMag)
    g90[i] = 1.645*gSTD*1000

# plot
plt.figure()
plt.scatter(pMag,g90,label=r"|$\vec p$|")
plt.scatter(rMag,g90,label=r"|$\vec r$|")
plt.title("Part 2")
plt.legend()
plt.ylabel("Error on Ground (m)")
plt.xlabel("Distance (km)")
print("Plot")

# Answer to Part 3
print("Part 3")

T = [T0[i] for i in range(3)]

for i in range(32):
    # Step through values of T from 35.3 to 45.3
    T[0] = T0[0] + 10/31*i

    # SIMULATE
    g,gMag,Tg = geoMonteSim(r0, v0, T, JD, err, n)

    # Recalc p vector for graphing
    p = [r0[i] - Tg[i] for i in range(3)]
    pMag[i] = norm(p)

    # Take standard deviation, use that to find the 90% confidence interval in meters
    gSTD = std(gMag)
    g90[i] = 1.645*gSTD*1000

# plot
plt.figure()
plt.scatter(pMag,g90,label=r"|$\vec p$|")
plt.xticks(range(1300,2000,100))
plt.title("Part 3")
plt.legend()
plt.ylabel("Error on Ground (m)")
plt.xlabel("Distance (km)")
print("Plot")

# Answer to part 4
print("Part 4")

# Initialize the dictionary that will be passed into function with one error at a time
singleErr = {
    'time': 0,
    'lati': 0,
    'long': 0,
    'elev': 0,
    'loca': 0,
    'V_xy': 0,
    'V__z': 0,
    'sPos': 0,
    'sAng': 0,
}

# initialize empty dictionary to store individual g vectors
g90 = {}

for key in err:
    if 'V' in key:
        # if the error is a velocity error it needs the time error in order to take effect
        singleErr['time'] = err['time']

    # put proper err in the dict
    singleErr[key] = err[key]

    # Run sim with only that err
    gMag = geoMonteSim(r0, v0, T0, JD, singleErr, n)[1]

    # Find 90% percent confidence interval, store in g90 dict
    gSTD = std(gMag)
    g90[key] = 1.645*gSTD*1000

    # Restore singleErr dict back to all zeros
    singleErr[key] = 0
    if 'V' in key:
        singleErr['time'] = 0

# Pull out names and values of 90% confidence intervael for graphing
names = list(g90.keys())
values = list(g90.values())

# Graph
plt.figure()
plt.bar(names, values)
plt.title("Part 4")
plt.ylabel("Error on Ground (m)")
plt.xlabel("Error Source")
print("Plot")

print("Part 5")
g90RMS = norm(values)
print(g90RMS)
print("This is usually higher than the one calculated in part 4.")
print("This is due to me passing in the time error along with the")
print("velocity error to actually make the velocity error apply.")

print("Part 6")

# calculate radius of circle fov presents on the ground assuming earth's curvature negligble
fov = 1e-2 # deg
fovRad = 1000000*m.tan(m.radians(1e-2)) # m

# Subtract the 90% confidence interval of the error on the ground since it counteracts the
# size of the circle you can be confident you are seeing
radius = fovRad - gRad

# Find the area of the square who's diagonal is the radius of the circle
# since the side length would be 2*sqrt(2)*R/2 the area is (sqrt(2)*R)^2
# which is equal to 2*R^2
squareArea = 2*radius**2

print("The area of the square we can be 90% confident our target is in is:",squareArea,"m^2")

plt.show()