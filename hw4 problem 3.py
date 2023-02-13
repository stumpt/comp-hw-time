import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.backends.backend_pdf import PdfPages
from numpy.linalg import inv
from functions import *

rcParams.update({'figure.autolayout': True})

#this simply suppresses scientific notation, just useful for debugging
np.set_printoptions(suppress=True)

#variable definitions
#unit set: kip, in, in2, kip*in, ksi

#define number of elements, with number of nodes 1 more than that
numElems = 3
equalLength = True #input whether elements will be equal length
numNodes = numElems + 1

#E (ksi), A (in2), L (in), P (kip)
E = 30000 #30E3
L = 36 #3 ft
P = -1.8 #1800lbf -> 1.8 kip, down
mConstant = (1.5)/(-36) #from geometry of cone
bConstant = 2
pi = np.pi

#element length = total length/number of elements if equal length elements
if equalLength == True:
    elementL = L/numElems

print("element length =" + str(elementL))

#displacements in mm
uList = [None]*numNodes #list of unknowns (node displacements)
#we can input our known values
uList[0] = 0 #fixed end at node 1
print("u list: ")
print(uList)

#forces
forceList = [0]*numNodes
forceList[0] = None
forceList[len(forceList) - 1] = P
print("force list: ")
print(forceList)


#empty ft matrix defined as 2 rows, (number of elements) columns
#element i's ft matrix is within column (i-1) -> element 1 ft values in column 0.
kMatrix = [0]*numElems
print("empty k matrix:")
print(kMatrix)

#begin by creating an array of zeroes with dimensions of [nodes x nodes]
gStiff = np.zeros((numNodes, numNodes))

#loop through each element
for i in range(numElems):
    print("i = ")
    print(i)
    """
    (x+(i*elementL))**2) -> **2 is squared (from surface traction function)
    x+ i*elementL adjusts local element x axis to global x
    {1-x/L, x/L} corresponds to [N]T
    """
    kStar = lambda x: ((mConstant*(x+(i*elementL)) + bConstant)**2) * ((E*pi)/(elementL**2))

    # integrate tstar from 0 to element length
    kVar, error = integrate.quad(kStar, 0, elementL)
    
    print("kVar = ")
    print(kVar)

    #send values to global ft matrix, allowing for easy access later
    #Ft(i) is accessable by at {ft[0,i], ft[1,i]}
    kMatrix[i] = kVar

    #each loop add bar constant to existing value in gstiff.
    gStiff[i,i] += kVar
    gStiff[i,(i+1)] += -kVar
    gStiff[(i+1),i] += -kVar
    gStiff[(i+1),(i+1)] += kVar
    
print("\nk Matrix = ")
print(kMatrix)

print("\ndisplacements = ")
print(uList)

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
    #first case: displacement BC = 0. In this case, set combined force values to 0, and alter 
    #gstiff to reflect the changes.
    if uList[j] == 0:
        print("u = 0 at row: " + str(j))
        forceList[j] = 0 
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
        forceList[j] = uList[j]
        for c in range(numNodes):
            if (gStiff[j,c] > 0) == True:
                gStiff[j,c] = 1
            else:
                gStiff[j,c] = 0

print("gstiff = ")
print(gStiff)

print("\nmod forces matrix: ")
print(forceList)

#multiply inverse global stiffness matrix by force matrix to get displacement
gStiffinv = inv(gStiff)
u = np.matmul(gStiffinv,forceList)

#conver to floats for easier manipulation
uFinal = u.astype(float)

#find stresses

stress = np.zeros([numElems, 1])

for z in range(numElems):
    stress[z] = (E/elementL)*(u[(z+1)]-u[z])

print(stress)
print(gStiffinv)

np.set_printoptions(
   formatter={'float':"{:0.8f}".format}, linewidth=130)

printNodalDisplacements(uFinal, "in")

#this line of code sets decimal places for final u.
np.set_printoptions(suppress=True,
   formatter={'float':"{:0.6f}".format}, linewidth=130)

printStresses(stress, "ksi")

#make output pdf
pdf = PdfPages('Hw4_problem3_outputs.pdf')

#make figure for displacement field
dispFig = plt.figure("Displacement Field")

for i in range(numElems):
    xVals = np.linspace((i*elementL), ((i+1)*(elementL)), 100)
    def y(x): return ((u[i+1] - u[i])/(elementL))*(x-(i*elementL)) + u[i]
    yVals = list(map(y, xVals))
    plt.plot(xVals, yVals)
    #plot labels
    plt.title("Problem 3: Displacement Field")
    plt.ylim(np.min(u)*2, 0)
    plt.xlim(0, L)
    plt.ylabel("u [in]")
    plt.xlabel("x [in]")

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
    plt.title("Problem 3: Stress Field")
    plt.ylim(np.min(stress)*2, 0)
    plt.xlim(0, L)
    plt.ylabel("\u03C3 [ksi]")
    plt.xlabel("x [in]")

pdf.savefig(stressFig)

#show plots
plt.show()
pdf.close()
