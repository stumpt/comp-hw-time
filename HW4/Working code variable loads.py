import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from numpy.linalg import inv
from functions import *

#this simply suppresses scientific notation, just useful for debugging
np.set_printoptions(suppress=True)

#variable definitions
#unit set: kip, in, in2, kip*in, ksi

#define number of elements, with number of nodes 1 more than that
numElems = 3
equalLength = True #input whether elements will be equal length
numNodes = numElems + 1

#E (MPa), A (mm2), L (mm), t*a (N/mm)
E = 70000 #70GPa
A = 300
L = 900 #5 ft
tb = 4.5 #600 lbf/in = 0.6 kip/in
mConstant = -0.005 #y=mx+b -> b= 4.5 (0,4.5) -> m= -4.5/900 = -0.005 (900,0)
bConstant = 4.5

#element length = total length/number of elements if equal length elements
if equalLength == True:
    elementL = L/numElems

print("element length =" + str(elementL))

barConstant = (E*A)/elementL

print("bar constant = " + str(barConstant))

#displacements in mm
uList = [None]*numNodes #list of unknowns (node displacements)
#we can input our known values
uList[0] = 0 #fixed end at node 1
uList[len(uList) - 1] = 0.02 
print("u list:")
print(uList)

#forces
forceList = [0]*numNodes
forceList[0] = None
forceList[len(forceList) - 1] = None
print("force list: ")
print(forceList)


#empty ft matrix defined as 2 rows, (number of elements) columns
#element i's ft matrix is within column (i-1) -> element 1 ft values in column 0.
ft = np.zeros([2, numElems])
print("empty ft:")
print(ft)

#begin by creating an array of zeroes with dimensions of [nodes x nodes]
gStiff = np.zeros((numNodes, numNodes))

#loop through each element
for i in range(numElems):
    #print("i = ")
    #print(i)
    """
    (x+(i*elementL))*m + b) -> surface traction function
    x+ i*elementL adjusts local element x axis to global x
    {1-x/L, x/L} corresponds to [N]T
    """
    tStar1 = lambda x: (mConstant*(x+(i*elementL)) + bConstant) * (1 - (x/elementL))
    tStar2 = lambda x: (mConstant*(x+(i*elementL)) + bConstant) * (x/elementL)

    # integrate tstar from 0 to element length
    ft1, error = integrate.quad(tStar1, 0, elementL) 
    ft2, error = integrate.quad(tStar2, 0, elementL)
    
    print("ft1 & 2 for element " + str(i))
    print(ft1)
    print(ft2)

    #send values to global ft matrix, allowing for easy access later
    #Ft(i) is accessable by at {ft[0,i], ft[1,i]}
    ft[[0],[i]] = ft1
    ft[[1],[i]] = ft2

    #each loop add bar constant to existing value in gstiff.
    gStiff[i,i] += barConstant
    gStiff[i,(i+1)] += -barConstant
    gStiff[(i+1),i] += -barConstant
    gStiff[(i+1),(i+1)] += barConstant
    
print("\nft = ")
print(ft)
ftFlat = ft.flatten('F')
print("\nflattened Ft: ")
print(ftFlat)
ftNumElems = len(ftFlat)
print("\nlength flat ft: ")
print(ftNumElems)

globalft = np.zeros((numNodes,1))
# this loop combines ft columns into one matrix
for k in range(numNodes):
    if k == 0:
        globalft[k] = ft[0,0]
    elif k == (numNodes - 1):
        globalft[k] = ftFlat[ftNumElems-1]
    else:
        globalft[k] = ftFlat[2*k] + ftFlat[(2*k) - 1]
        
print("\nGlobal ft matrix: ")
print(globalft)

totalForces = np.zeros_like(globalft)
for i1 in range(numNodes):
    if forceList[i1] != None:
        totalForces[i1] = forceList[i1] + globalft[i1]
    else:
        totalForces[i1] = globalft[i1]

print("\ntotal forces matrix: ")
print(totalForces)
print("\ndisplacements = ")
print(uList)
print("\nConstant [kip/in] = ")
print(barConstant)

check_Symmetric(gStiff)
print(gStiff)

"""
This code is "method 3," begins by checking if a displacement is 0. if it is, it sets the 
corresponding force value to 0. Then, it checks for positivity/negativity within the row
that the 0 displacement corresponds to. Positive values get set to 1, while negatives
get set to 0. It also catches None (unknown) values by just incrementing the loop value.
For any other known u values, it sets the force equal to that value and proceeds to alter
the stiffness matrix in the same fashion.
"""
for j in range(numNodes):
    print(str(j))
    #first case: displacement BC = 0. In this case, set combined force values to 0, and alter 
    #gstiff to reflect the changes.
    if uList[j] == 0:
        print("u = 0 at row: " + str(j))
        totalForces[j] = 0 
        for c in range(numNodes):
            if (gStiff[j,c] > 0) == True:
                gStiff[j,c] = 1
            else:
                gStiff[j,c] = 0
    elif uList[j] == None:
        print("none detected at row " + str(j))
    else: #when we have a known, non-zero displacement boundary condition
        #alter global stiffness to all zeros and one 1, and set force equal to displacement BC
        print("row " + str(j) + " altered")
        print(uList[j])
        totalForces[j] = uList[j]
        for c in range(numNodes):
            if (gStiff[j,c] > 0) == True:
                gStiff[j,c] = 1
            else:
                gStiff[j,c] = 0

print("gstiff = ")
print(gStiff)

print("\ntotal forces matrix: ")
print(totalForces)

#multiply inverse global stiffness matrix by force matrix to get displacement
gStiffinv = inv(gStiff)
print("\n inverted global stiffness: ")
print(gStiffinv)
u = np.matmul(gStiffinv,totalForces)

#conver to floats for easier manipulation
uFinal = u.astype(float)

#find stresses

stress = np.zeros([numElems, 1])

for z in range(numElems):
    stress[z] = (E/elementL)*(u[(z+1)]-u[z])

print(stress)



#this line of code sets 4 decimal places for final u.
np.set_printoptions(suppress=True,
   formatter={'float':"{:0.6f}".format}, linewidth=130)

printNodalDisplacements(uFinal, "mm")
printStresses(stress, "MPa")

#make output pdf
pdf = PdfPages('Hw4_problem2_outputs.pdf')

#make figure for displacement field
dispFig = plt.figure("Displacement Field")

for i in range(numElems):
    xVals = np.linspace((i*elementL), ((i+1)*(elementL)), 100)
    def y(x): return ((u[i+1] - u[i])/(elementL))*(x-(i*elementL)) + u[i]
    yVals = list(map(y, xVals))
    plt.plot(xVals, yVals)
    #plot labels
    plt.title("Problem 2 Part 1: Displacement Field")
    plt.ylim(0, (np.max(u))*1.1)
    plt.xlim(0, L)
    plt.axhline(0, color="black", linestyle = "--")
    plt.axvline(0, color="black", linestyle = "--")
    plt.ylabel("u [mm]")
    plt.xlabel("x [mm]")

pdf.savefig(dispFig)

#make figure for stress field
stressFig = plt.figure("Stress Field")

for i in range(numElems):
    #x values sets each interval's values
    xValues = np.linspace((i*elementL), ((i+1)*(elementL)), 100)
    def f(x): return stress[i]
    yValues = list(map(f, xValues))
    plt.plot(xValues, yValues)
    #plot labels
    plt.title("Problem 2 Part 1: Stress Field")
    plt.ylim((np.min(stress)*1.3), (np.max(stress))*1.1)
    plt.xlim(0, L)
    plt.axhline(0, color="black", linestyle = "--")
    plt.axvline(0, color="black", linestyle = "--")
    plt.ylabel("\u03C3 [MPa]")
    plt.xlabel("x [mm]")

pdf.savefig(stressFig)

#show plots
plt.show()
pdf.close()
