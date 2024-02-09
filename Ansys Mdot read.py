# Reads all of the mass flow rate report files given by ANSYS for the CARILLON liquid rocket engine
# Written by Joel Surfleet 2/6/24 for CPSS

import os
import glob
import pandas as pd

os.chdir(r"C:\Users\joels\ox_sim_feb5_files\dp0\FFF-1\Fluent")  # Paste in the directory where report files are stored (should be called 'Fluent')

print("Loaded into:",os.getcwd())   # Confirms proper file path
print()

nameList = glob.glob("*.out")   # Gets a list of all of the report files

print(nameList)                 # Prints those names so you can be sure you don't have anything extra along for the ride
print()

size = len(nameList)    # Finds the number of data points (should be 20 since we have 20 like doublets)

dataList = [None] * size

for name in nameList:   # This loop opens all of the .out files, reads the final value and writes it to then new .txt file

    file = open(name,'r')
    data = file.read()
    file.close()
    
    data = data.splitlines()

    data = (float(data[-1][4:]) + float(data[-1][4:])) / 2      # Collects the cumulative average of the last 2 data points

    index = name[10:12].split("-")                              # Figures out which hole is being observed

    index = int(index[0]) - 1                                   # Makes that hold number into an index

    dataList[index] = data                                      # passes Mdot values into an array in the proper order

total = 0   # initialize total for averaging

maxi = 0
maxM = 0
mini = 0
minM = 10

for i in range(size):   
    total += dataList[i]            # Sums the total mass flow rate

    if abs(dataList[i]) > abs(maxM):          # Stores the max value and index
        maxM = dataList[i]
        maxi = i

    if abs(dataList[i]) < abs(minM):          # Stores the min value and index
        minM = dataList[i]
        mini = i

avgMdot = total / size  # Finds the average value of the Mdot accross the holes

percentDiff = [None] * size

for i in range(size):
    percentDiff[i] = ((dataList[i] - avgMdot) / avgMdot) * 100     # Caclulates percent different from average

data = {
    "Hole Number": range(1,size + 1),
    "Mass Flow Rate": dataList,
    "Percent Different from Average": percentDiff
}

df = pd.DataFrame(data)

df.to_excel('Mass Flow Rate.xlsx',index=False)