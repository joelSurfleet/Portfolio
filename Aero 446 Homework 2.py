# Aero 446 HW2

print("\nProblem 2")

def MGA(sat,states):

    # Implementation of the SAWE A-3 mass properties standard
    # Used to estimate the mass growth allowance of sattelite subsystems
    # Make sure you have my SAWE_A3.txt file and copy the path below

    # Inputs
    #   sat: a dictionary containing current mass values of the subsystems
    #   states: a dictionary containg development stages of the subsystems
    # make sure you use the naming conventions found in the first row and column of the txt

    # Outputs
    # Adds new variables to the sat dictionary 
    # calls them "subName"_min and "subName"_max for each subsystem

    f = open(r"C:\Users\joels\Documents\Python\SAWE_A3.txt")
    SAWE = f.read()
    f.close()

    SAWE = SAWE.splitlines()
    SAWE = [SAWE[i].split() for i in range(6)]

    for sub in states:

        for j in range(len(SAWE[0])):
            if SAWE[0][j] == sub:
                id = j
                break

        if states[sub] == "NTE":
            sat[sub+"_min"] = sat[sub]
            sat[sub+"_max"] = sat[sub]

        else:
            for i in range(1,6):
                if SAWE[i][0] == states[sub]:
                    growth = SAWE[i][j].split('-')
                    sat[sub+"_min"] = sat[sub] * (1 + (float(growth[0]) / 100))
                    sat[sub+"_max"] = sat[sub] * (1 + (float(growth[1]) / 100))

sattelite = {
    # Mass units in kg
    "ElecS" : 350,
    "ElecM" : 0,
    "ElecL" : 0,
    "Struc" : 270,
    "Therm" : 170,
    "Props" : 200,
    "Batts" : 125,
    "Solar" : 180,
    "Wires" : 75,
    "Mechs" : 10
}

states = {
    "ElecS" : "NTE",
    "ElecM" : "NTE",
    "ElecL" : "NTE",
    "Struc" : "Layout",
    "Therm" : "Estimated",
    "Props" : "Layout",
    "Batts" : "Estimated",
    "Solar" : "Estimated",
    "Wires" : "Estimated",
    "Mechs" : "Existing"
}

MGA(sattelite,states)

MGA_min = 0
MGA_max = 0
massInit= 0

for sub in sattelite:
    if '_min' in sub:
        MGA_min = MGA_min + sattelite[sub]
    elif '_max' in sub:
        MGA_max = MGA_max + sattelite[sub]
    else:
        massInit= massInit+ sattelite[sub]

print("Min Mass Growth Allowed for whole system:")
print(MGA_min,"kg")
print("Max Mass Growth Allowed for whole system:")
print(MGA_max,"kg")
print()

print("Problem 3")
print("Percent of total mass for growth, Minimum:")
print(MGA_min / (massInit + MGA_min),"%")
print("Percent of total mass for growth, Minimum:")
print(MGA_max / (massInit + MGA_max),"%")
print()