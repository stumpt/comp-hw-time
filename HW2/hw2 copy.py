import numpy as np

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
            print ("The matrix is symmetric.")
        else:
            print ("Error 1: The matrix is not symmetric.")
    else:
        print ("Error 2: The matrix is not symmetric.")

#declare arrays for all givens, which allows for looping, and record lengths
k = np.matrix([[k1], [k2], [k3]])
kLength = len(k)

P = np.matrix([[P1], [P2], [P3], [P4]])
PLength = len(P)

u = np.matrix([[u1], [u2], [u3], [u4]])
uLength = len(u)

#define iterator for loops
i=0
i1=0
i2=0

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


while j < uLength:
    if u[j] == 0:
        P[j] = 0 
        while k < numNodes:
            if (gStiff[j,k] > 0) == True:
                gStiff[j,k] = 1
            else:
                gStiff[j,k] = 0
            k+=1




#gStiff[0,1] = k1

print(gStiff)
print (P)
