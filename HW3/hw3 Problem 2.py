import numpy as np
from numpy.linalg import inv

#this simply suppresses scientific notation, just useful for debugging
np.set_printoptions(suppress=True)

#variable definitions
#units: N, mm, MPa

#E (MPa), A (mm^2), and L (mm)
E = 200000
A1 = 2000
A2 = 500
L1 = 600
L2 = 500

#concentrated/unknown loads in N - None indicates unknown
P1 = None
P2 = 200000
P3 = None


#displacements in mm
u1 = 0
u2 = None
u3 = 0.15


#define number of elements, with number of nodes 1 more than that
numElems = 2
numNodes = numElems + 1

#define function to check if an array/matrix is symmetric
def check_Symmetric(array):
    transposed = array.transpose()
    if array.shape == transposed.shape:
        if (array == transposed).all():
            print ("The matrix is symmetric. \n")
        else:
            print ("Error 1: The matrix is not symmetric. \n")
    else:
        print ("Error 2: The matrix is not symmetric. \n")

#declare arrays for all givens, which allows for looping, and record lengths
areas = np.matrix([[A1], [A2]])
ALength = len(areas)

lengths = np.matrix([[L1], [L2]])
lLength = len(lengths)

forces = np.matrix([[P1], [P2], [P3]])
PLength = len(forces)

displacements = np.matrix([[u1], [u2], [u3]])
uLength = len(displacements)

stress = np.zeros((numElems, 1))

barConstants = np.matrix([[(E*A1)/L1], [(E*A2/L2)]])

""" Debug Code
print("Areas = ")
print(areas)
print("Lengths = ")
print(lengths)
print("displacements = ")
print(displacements)
print("Bar constants = ")
print(barConstants)"""

#define iterator for loops
i=0
j=0
c=0
z=0

#create global stiffness

#begin by creating an array of zeroes with dimensions of [nodes x nodes]
gStiff = np.zeros((numNodes, numNodes))

#loop through each element (numElems times), each time assigning element values to +- k[i] values
#+= ensures new values are added onto previous instead of being rewritten
while i < numElems:
    gStiff[i,i] += barConstants[i]
    gStiff[i,(i+1)] += -barConstants[i]
    gStiff[(i+1),i] += -barConstants[i]
    gStiff[(i+1),(i+1)] += barConstants[i]
    #iterate counter
    i+=1

check_Symmetric(gStiff)
#print(gStiff)

"""
This code is "method 3," begins by checking if a displacement is 0. if it is, it sets the 
corresponding force value to 0. Then, it checks for positivity/negativity within the row
that the 0 displacement corresponds to. Positive values get set to 1, while negatives
get set to 0. It also catches None (unknown) values by just incrementing the loop value.
For any other known u values, it sets the force equal to that value and proceeds to alter
the stiffness matrix in the same fashion.
"""
while j < uLength:
    if displacements[j] == 0:
        c=0
        print("u = 0 at row: " + str(j))
        forces[j] = 0 
        while c < numNodes:
            if (gStiff[j,c] > 0) == True:
                gStiff[j,c] = 1
            else:
                gStiff[j,c] = 0
            c+=1
            print(c)
        j+=1
    elif displacements[j] == None:
        print("none detected at row " + str(j))
        j+=1
    else:
        c=0
        print("row " + str(j) + " altered")
        print(displacements[j])
        forces[j] = displacements[j]
        while c < numNodes:
            if (gStiff[j,c] > 0) == True:
                gStiff[j,c] = 1
            else:
                gStiff[j,c] = 0
            c+=1
            print(c)
        j+=1

#multiply inverse global stiffness matrix by force matrix to get displacement
gStiffinv = inv(gStiff)
u = np.matmul(gStiffinv,forces)

#conver to floats for easier manipulation
uFinal = u.astype(float)

#find forces

while z < numElems:
    stress[z] = (E/lengths[z])*(u[(z+1)]-u[z])
    z+=1       

#print(gStiff)
#print(gStiffinv)


#this line of code sets 4 decimal places for final u.
np.set_printoptions(suppress=True,
   formatter={'float':"{:0.4f}".format}, linewidth=130)

print("\nThe nodal displacements, in [mm], ordered from lowest to highest numbered node are:\n")
print(uFinal)
print("\n")

print("The element axial stresses, in [MPa], ordered from lowest to highest numbered element are:\n ")
print(stress)
print("\nNote: a (+) value indicates the element is in tension, while a (-) value in compression.\n")

