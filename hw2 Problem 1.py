import numpy as np
from numpy.linalg import inv

#this simply suppresses scientific notation, just useful for debugging
np.set_printoptions(suppress=True)

#variable definitions

#spring constants in lb/in
k1 = 5000 
k2 = 5000
k3 = 5000

#concentrated/unknown loads in lbf
P1 = None
P2 = -1000
P3 = 0 
P4 = 4000

#displacements in inches
u1 = 0
u2 = None
u3 = None
u4 = None


#define number of elements, with number of nodes 1 more than that
numElems = 3
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
k = np.matrix([[k1], [k2], [k3]])
kLength = len(k)

P = np.matrix([[P1], [P2], [P3], [P4]])
PLength = len(P)

u = np.matrix([[u1], [u2], [u3], [u4]])
uLength = len(u)

PElems = np.zeros((numElems, 1))

#define iterator for loops
i=0
j=0
c=0
z=0

#create global stiffness

#begin by creating an array of zeroes with dimensions of [nodes x nodes]
gStiff = np.zeros((numNodes, numNodes))

#loop through each element (kLength times), each time assigning element values to +- k[i] values
#+= ensures new values are added onto previous instead of being rewritten
while i < kLength:
    gStiff[i,i] += k[i]
    gStiff[i,(i+1)] += -k[i]
    gStiff[(i+1),i] += -k[i]
    gStiff[(i+1),(i+1)] += k[i]
    #iterate counter
    i+=1

check_Symmetric(gStiff)

"""
This code is "method 3," begins by checking if a displacement is 0. if it is, it sets the 
corresponding force value to 0. Then, it checks for positivity/negativity within the row
that the 0 displacement corresponds to. Positive values get set to 1, while negatives
get set to 0.
"""
while j < uLength:
    if u[j] == 0:
        c=0
        P[j] = 0 
        while c < numNodes:
            if (gStiff[j,c] > 0) == True:
                gStiff[j,c] = 1
            else:
                gStiff[j,c] = 0
            c+=1
    j+=1



#multiply inverse global stiffness matrix by force matrix to get displacement
gStiffinv = inv(gStiff)
u = np.matmul(gStiffinv,P)

#conver to floats for easier manipulation
uFinal = u.astype(float)

#find forces
"""The absolute value of the force within each element is the same, thus that value can 
be found by taking the absolute value of kx*ux - kx * u(x+1),where x is any element. 
If the non absolute value yielded from that equation is negative, the 
element is in tension. Likewise, if it is positive, the element is in compression.
This value will be multiplied by -1 to match standard sign conventions."""

while z < numElems:
    PElems[z] = (k[z]*uFinal[z] - k[z]*uFinal[(z+1)])*-1
    z+=1         

#this line of code sets 4 decimal places for final u.
np.set_printoptions(suppress=True,
   formatter={'float':"{:0.4f}".format}, linewidth=130)

print("The nodal displacements, in [in], ordered from lowest to highest numbered node are:\n")
print(uFinal)
print("\n")

print("The forces in the springs, in [lbf], ordered from lowest to highest numbered element are:\n ")
print(PElems)
print("\nNote: a (+) value indicates the element is in tension, while a (-) value in compression.\n")
