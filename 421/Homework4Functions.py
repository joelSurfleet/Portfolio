# Homework4Functions

import math as m
from numpy.random import standard_normal

def lst(lon, JD):
    # CT2LST  Find the local sidereal time from the local Civilian Time in
    # degrees and Julian Date
    
    # lst = ct2lst(lon, JD)
    # Returns the local sidereal time in degrees, lst, given the local civilian
    # time in degrees, lon, and the Julian Date, JD.  The local Civilian Time
    # is the same as the East Longitude of the current location of interest on
    # earth.
    
    # From Curtis
    # Author: Eric A. Mehiel
    # Date: December 19th, 2006

    # Translated to Python by Joel

    J0 = m.floor(JD + 0.5) - 0.5
    UT = (JD - J0)*24.0
    T0 = (J0 - 2451545.0)/36525.0
    c = [100.4606184, 36000.77004, 0.000387933, -2.583e-8, 360.98564724]
    g0 = c[0] + c[1]*T0 + c[2]*T0**2 + c[3]*T0**3

    g0 = g0 - 360.0*m.floor(g0/360.0)

    gst = g0 + c[4]*UT/24.0
    lst = gst + lon

    lst = lst - 360.0*m.floor(lst/360.0)

    return lst

def lla2eci(pos, JD):
    # Converts a lat long alt coordinate into ECI
    # Take from Dr. Mehiel
    # Originally written in MatLab, translated to python

    lat = pos[0]
    lon = pos[1]
    alt = pos[2]

    Re = 6378.137
    gst = lst(lon, JD) #deg
    theta = (gst) * m.pi/180.0 #rad

    r = (Re + alt)*m.cos(m.radians(lat)) # km

    pos = [r*m.cos(theta), r*m.sin(theta), (alt+Re)*m.sin(m.radians(lat))] # km

    return pos

def norm(a,n=2):
    # Takes the norm of a vector, default is 2-norm if none is specified

    if n == 'inf':
        aMag = 0
        for val in a:
            if val > aMag:
                aMag = val

    else:
        sum = 0
        for val in a:
            sum = sum + val**n
        aMag = sum ** (1/n)

    return aMag

def dot(a,b):
    # Dot product of two vectors

    tot = 0

    if len(a) != len(b):
        print("Vectors Must be the same length")
        return 0
    
    for i in range(len(a)):
        tot = tot + a[i]*b[i]
    
    return tot

# Store the square value of earth's radius for speed of execution
Re2 = 6378.137**2

def earthIntersect(pHat,r):
    # Finds the point that a vector intersects the earth's surface
    # Outputs a Vector in ECI from the center of earth, and the distance from
    # r to the tip of that vector.

    rMag = norm(r)
    
    d1 = -dot(r,pHat) + (dot(r,pHat)**(2) - (rMag**2 - Re2)) ** (1/2)
    d2 = -dot(r,pHat) - (dot(r,pHat)**(2) - (rMag**2 - Re2)) ** (1/2)

    if abs(d1) > abs(d2):
        p = [d1*p for p in pHat]
        T0 = [r[i]+p[i] for i in range(len(r))]
        return T0, d1
    else:
        p = [d2*p for p in pHat]
        T0 = [r[i]+p[i] for i in range(len(r))]
        return T0, d2
    
def pointingError(vector,std):
    # Finds the pointing error of a vector given a standard deviation
    # Applies to rotation errors in x, y and z

    xErr = std*standard_normal()
    yErr = std*standard_normal()
    zErr = std*standard_normal()

    Cx = basicRotation('x',xErr)
    Cy = basicRotation('y',yErr)
    Cz = basicRotation('z',zErr)

    vector = [sum(Cx[0])*vector[0], sum(Cx[1])*vector[1], sum(Cx[2])*vector[2]]
    vector = [sum(Cy[0])*vector[0], sum(Cy[1])*vector[1], sum(Cy[2])*vector[2]]
    vector = [sum(Cz[0])*vector[0], sum(Cz[1])*vector[1], sum(Cz[2])*vector[2]]

    return vector

def basicRotation(axis,deg):
    # Basic rotation matrices in one axis at a time
    # inputs: axis - string argument either 'x','y', or 'z'
    # inputs: deg - rotation in degrees
    # outputs: direction cosine matrix

    theta = m.radians(deg)

    if axis == 'x':
        C = [[1,0,0], [0,m.cos(theta),m.sin(theta)], [0,-m.sin(theta),m.cos(theta)]]
    elif axis == 'y':
        C = [[m.cos(theta),0,-m.sin(theta)], [0,1,0], [m.sin(theta),0,m.cos(theta)]]
    elif axis == 'z':
        C = [[m.cos(theta),m.sin(theta),0], [-m.sin(theta),m.cos(theta),0], [0,0,1]]

    return C

def geoMonteSim(R0, V0, T0, JD, err, n):
    # Runs a monte-carlo simulation to simulate pointing error of a sattelite.

    # Find P, then normalize it to get pHat
    p = [R0[i] - T0[i] for i in range(3)]
    pMag = norm(p)
    pHat = [p[i] / pMag for i in range(3)]
    
    # Initialize vectors that will be calculated and stored n number of times
    rErr = [0,0,0]
    tErr = [0,0,0]
    pErr = [0,0,0]
    g    = [[0,0,0] for i in range(n)]
    gMag = [[0,0,0] for i in range(n)]
    Tg = [0,0,0]

    for i in range(n):

        # Calcualte position with random positional error
        rErr[0] = R0[0] + err['loca']*standard_normal() + err['time']*(err['V_xy'])*standard_normal() + err['sPos']*standard_normal()
        rErr[1] = R0[1] + err['loca']*standard_normal() + err['time']*(err['V_xy'])*standard_normal() + err['sPos']*standard_normal()
        rErr[2] = R0[2] + err['loca']*standard_normal() + err['time']*(err['V__z'])*standard_normal()

        # Calculate target position with error in lat-long-alt
        tErr[0] = T0[0] + err['lati']*standard_normal()
        tErr[1] = T0[1] + err['long']*standard_normal()
        tErr[2] = T0[2] + err['elev']*standard_normal()

        # Calculate time error
        timeErr = (err['time']/86400)*standard_normal()
        lstErr = lst(tErr[1], JD + timeErr)
        JDloc = JD + lstErr/360

        # convert target position into ECI
        tErr = lla2eci(tErr, JDloc)

        # Sum total target position for averaging later
        Tg = [Tg[j] + tErr[j] for j in range(3)]

        # Calculate error in p vector
        pErr = pointingError(pHat,err['sAng'])
        pErr = [pErr[j] * pMag for j in range(3)]

        # Aggregate errors into an error on the ground
        g[i] = [pErr[j]+rErr[j]-tErr[j] for j in range(3)]
        gMag[i] = norm(g[i])

    # Average value of target vector
    Tg = [Tg[i]/n for i in range(3)]

    return g,gMag,Tg

def std(list):
    # Calculates standard deviation of a list

    n = len(list)
    mean = sum(list)/n

    sumDiffSquares = 0

    for val in list:
        sumDiffSquares = sumDiffSquares + ((val - mean) ** 2)

    return ((sumDiffSquares/(n-1)) ** (1/2))