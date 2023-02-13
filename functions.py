import numpy as np
from numpy.linalg import inv

#function to check if an array/matrix is symmetric
def check_Symmetric(array):
    transposed = array.transpose()
    if array.shape == transposed.shape:
        if (array == transposed).all():
            print ("The matrix is symmetric. \n")
        else:
            print ("Error 1: The matrix is not symmetric. \n")
    else:
        print ("Error 2: The matrix is not symmetric. \n")

#function to create [N]
def create_N(X, L):
    N1 = (1/L**3) * ((2)*(X**3) - (3)*(X**2)*L + L**3)
    N2 = (1/L**3) * ((X**3)*L - (2)*(X**2)*(L**2) + X*(L**3))
    N3 = (1/L**3) * ((-2)*(X**3) + (3)*(X**2)*L)
    N4 = (1/L**3) * ((X**3)*L - (X**2)*(L**2))
    Nmatrix = np.matrix([[N1], [N2], [N3], [N4]])
    NmatrixT = np.transpose(Nmatrix) #transposed -> [N]T

def printNodalDisplacements(displacementArray, unit: str):
    print('\u2588'*100)
    print("The nodal displacements, in [" + unit + "], ordered from lowest to highest numbered node are:\n")
    print(displacementArray)
    print('\u2588'*100)
    


def printGivenNodeDisplacement(node,displacementArray, unit: str):
    print('\u2588'*100)
    print("\nThe displacement of point B (node " + str(node) + ", in [" + unit + "], is given as:\n")
    print(str(displacementArray[node - 1]))
    print('\u2588'*100)

def printSpringForce(springForce, unit: str):
    print('\u2588'*100)
    print("\nThe force in the spring, in [" + unit + "], is:\n")
    print(springForce)
    print("\nNote: a (+) value indicates the element is in tension, while a (-) value in compression.\n")
    print('\u2588'*100)

def printStresses(stressMatrix, unit:str):
    print('\u2588'*100)
    print("\nThe element axial stresses, in [" + unit + "], ordered from lowest to highest numbered element are:\n ")
    print(stressMatrix)
    print("\nNote: a (+) value indicates the element is in tension, while a (-) value in compression.\n")
    print('\u2588'*100)

