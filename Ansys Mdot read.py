# Reads all of the mass flow rate report files given by ANSYS for the CARILLON liquid rocket engine
# Written by Joel Surfleet 2/16/24 for CPSS

# Takes report-def.out files from ansys mass flow rate surface tracking and outputs an excel file with relevant data
# Copy in file path for report files
# The program will output "Mass Flow Rate.xlsx" to the same directory as the python script
# or will OVERWRITE a file with that name. Can't be run if that file is open in excel.

import pandas as pd
import math as m

f = open(r"C:\Users\joels\Carillon Flow Simulation_files\dp0\FFF\Fluent\report-def-0-rfile.out",'r') # Opens my folder for the mass flow rate data (ctrl+shft+c your file to get path)

dataString = f.read()
dataString = dataString.splitlines()[-1]
dataString = dataString.split(' ')          # Turns the read data for the last iteration into an array of strings

f.close()

size = len(dataString) - 1  # Gets the number of data points (subtracts one since the first number is the iteration #)

# Initialize all of the lists you need based on the number of data points
holeNum = [None] * size
fuelMdot = [None] * size
fuelMdotUnordered = [None] * size
fuelPercent = [None] * size
oxMdot = [None] * size
oxMdotUnordered = [None] * size
oxPercent = [None] * size
OFratio = [None] * size
OFpercent = [None] * size
resultant = [None] * size
resultantPercent = [None] * size

# This is the velocity out of the fuel and ox orifices, used for momentum calcs to find resultant angle
vf = 29.9
vO = 33.3

# Populate the hole number list correctly and turn string data into float data
for i in range(size):
    holeNum[i] = i
    fuelMdotUnordered[i] = abs(float(dataString[i+1]))

# Put the data in the correct order since ansys orders files really stupidly
fuelMdot[0:2]   = fuelMdotUnordered[0:2]
fuelMdot[2:10]  = fuelMdotUnordered[12:20]
fuelMdot[10:20] = fuelMdotUnordered[2:12]

# Get the datafile for the ox sim
f = open(r"C:\Users\joels\Carillon Flow Simulation_files\dp0\FFF-1\Fluent\report-def-0-rfile.out",'r')

dataString = f.read()
dataString = dataString.splitlines()[-1]
dataString = dataString.split(' ')

f.close()

# Turn it to float
for i in range(size):
    oxMdotUnordered[i] = abs(float(dataString[i+1]))

# Reorder
oxMdot[0:2]   = oxMdotUnordered[0:2]
oxMdot[2:10]  = oxMdotUnordered[12:20]
oxMdot[10:20] = oxMdotUnordered[2:12]

# Calculate the Mixture ratio and resultant angle after impingement
for i in range(size):
    resultant[i] = m.degrees(m.atan((oxMdot[i] * vO * m.sin(m.radians(60)) - fuelMdot[i] * vf * m.sin(0)) / (oxMdot[i] * vO * m.cos(m.radians(60)) + fuelMdot[i] * vf * m.cos(0))))
    OFratio[i] = oxMdot[i]/fuelMdot[i]

# Find the average for the lists (should be the same as the nominal value, but this takes into account the slight errors introduced by the simulation)
avgFuel = sum(fuelMdot)/len(fuelMdot)
avgOx = sum(oxMdot)/len(oxMdot)
avgOF = sum(OFratio) / size
avgResultant = sum(resultant) / size

# Make lists of the percent away from average for each thing we're tracking for each hole
for i in range(size):  
    fuelPercent[i] = (fuelMdot[i] - avgFuel) / avgFuel * 100
    oxPercent[i] = (oxMdot[i] - avgOx) / avgOx * 100
    OFpercent[i] = (OFratio[i] - avgOF) / avgOF * 100
    resultantPercent[i] = (OFratio[i] - avgOF) / avgOF * 100

# Put the data into a data frame with named columns
dataFrame = {
    "Hole Pair": range(size),
    "Fuel Mass Flow Rate": fuelMdot,
    "% from Average F": fuelPercent,
    "Ox Mass Flow Rate": oxMdot,
    "% Different from Average O": oxPercent,
    "O/F ratio": OFratio,
    "% from Average OF": OFpercent,
    "Resultant Angle": resultant,
    "% from Average d": resultantPercent,
}

# Output to excel in the same directory as where the PYTHON SCRIPT is located.
# Note that this can't be done if the program is open and it OVERWRITES ALL FORMATTING so be careful
df = pd.DataFrame(dataFrame)
df.to_excel('Mass Flow Rate.xlsx',index=False)