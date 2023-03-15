import numpy as np
from numpy.linalg import inv
import math as m

#function to check if an array/matrix is symmetric
def check_symmetric(array):
    transposed = array.transpose()
    if array.shape == transposed.shape:
        if (array == transposed).all():
            print ("The matrix is symmetric. \n")
        else:
            print ("Error 1: The matrix is not symmetric. \n")
    else:
        print ("Error 2: The matrix is not symmetric. \n")

#function to create [N]T
def create_NT(X, L):
    N1 = (1/L**3) * ((2)*(X**3) - (3)*(X**2)*L + L**3)
    N2 = (1/L**3) * ((X**3)*L - (2)*(X**2)*(L**2) + X*(L**3))
    N3 = (1/L**3) * ((-2)*(X**3) + (3)*(X**2)*L)
    N4 = (1/L**3) * ((X**3)*L - (X**2)*(L**2))
    return np.matrix([N1, N2, N3, N4])

#This simple function returns [k] for a given theta.
def create_truss_k(theta):
    return np.array([
    [(m.cos(theta))**2, m.cos(theta) * m.sin(theta), -((m.cos(theta))**2), -(m.cos(theta) * m.sin(theta))], 
    [m.cos(theta) * m.sin(theta), (m.sin(theta))**2, -(m.cos(theta) * m.sin(theta)), -((m.sin(theta))**2)], 
    [-((m.cos(theta))**2), -(m.cos(theta) * m.sin(theta)), (m.cos(theta))**2, m.cos(theta) * m.sin(theta)], 
    [-(m.cos(theta) * m.sin(theta)), -((m.sin(theta))**2), m.cos(theta) * m.sin(theta), (m.sin(theta))**2]
    ])

def create_beam_local_k(L):
    return np.array([
        [12, 6*L, -12, 6*L],
        [6*L, 4*(L**2), -6*L, 2*(L**2)],
        [-12, -6*L, 12, -6*L],
        [6*L, 2*(L**2), -6*L, 4*(L**2)]
    ])

def create_truss_stress_trig_matrix(theta):
    return np.array([
        [m.cos(theta), m.sin(theta), 0, 0],
        [0, 0, m.cos(theta), m.sin(theta)]
    ])

def get_element_length(x1, y1, x2, y2):
    x = (x2 - x1)**2
    y = (y2 - y1)**2
    return m.sqrt(x+y)


#This function takes in known boundary conditions and an array, translating inputs into
#usable arrays for the code.
def adjust_array(ArrayToAdjust, Given_BCs, tag: str):
    print(f"adjusting array to account for: {tag} boundary conditions")
	#first, count how many boundary conditions are given.
    loadRows = (np.shape(Given_BCs))[0]
    for i in range(loadRows):
    #adjusted node translates the given node to the corresponding array position for that node's forces.
        adjustedNode = int(2 * (Given_BCs[i,0] - 1))
        #recording the boundar conditions themselves is unnecessary, but facilitates easier code reading.
        nodeXCondition = Given_BCs[i,1]
        nodeYCondition = Given_BCs[i,2]
        #Check: does the given BC match the current array element? if it doesn't, then change the array to the given BC.
        #this is done like this such that it works for both forces and displacements
        if Given_BCs[i,1] != ArrayToAdjust[adjustedNode]:
            print(f"case 1 for BC {str(i + 1)}")
            ArrayToAdjust[adjustedNode] = nodeXCondition
        else:
            print(f"case 3 for BC {str(i + 1)}")

        if Given_BCs[i,2] != ArrayToAdjust[adjustedNode + 1]:
            print(f"Case 2 for BC {str(i + 1)}")
            ArrayToAdjust[adjustedNode + 1] = nodeYCondition
        else:
            print(f"case 4 for BC {str(i + 1)}")
    print(f"done adjusting {tag} array!")
    return ArrayToAdjust

    

#this function allows for a transplant point to be provided (new row & column), and changes the new array to match the old
#"block" within the old array. A block is essentially a 2x2 matrix of adjacent elements within a larger matrix.
def block_array_adjust(NewRow, NewColumn, NewArray, OldRow, OldColumn, OldArray,):
    NewArray[NewRow, NewColumn] += OldArray[OldRow, OldColumn]
    NewArray[NewRow, (NewColumn + 1)] += OldArray[OldRow, (OldColumn + 1)]
    NewArray[(NewRow + 1), NewColumn] += OldArray[(OldRow + 1),OldColumn]
    NewArray[(NewRow + 1), (NewColumn + 1)] += OldArray[(OldRow + 1), (OldColumn + 1)]

def check_even(num):
    return num % 2 == 0

def get_radius_from_I(I):
    return ((4*I)/m.pi)**0.25


    

def printNodalDisplacements(displacementArray, unit: str):
    print('\u2588'*100)
    print(f"The nodal displacements, in [{unit}], ordered from lowest to highest numbered node are:\n{displacementArray}")
    print('\u2588'*100)

    
def printTrussNodalDisplacements(displacementArray, unit: str):
    print('\u2588'*100)
    print(f"The nodal displacements (u,v), in [{unit}], ordered from lowest to highest numbered node are:\n{displacementArray}")
    print('\u2588'*100)
    
def printGivenNodeDisplacement(node,displacementArray, unit: str):
    print('\u2588'*100)
    print(f"\nThe displacement of point B (node {node}, in [{unit}], is given as:\n{displacementArray[node-1]}")
    print('\u2588'*100)

def printSpringForce(springForce, unit: str):
    print('\u2588'*100)
    print(f"\nThe force in the spring, in [{unit}], is:\n{springForce}")
    print("\nNote: a (+) value indicates the element is in tension, while a (-) value in compression.\n")
    print('\u2588'*100)

def printStresses(stressMatrix, unit:str):
    print('\u2588'*100)
    print(f"\nThe element axial stresses, in [{unit}], ordered from lowest to highest numbered element are:\n{stressMatrix}")
    print("\nNote: a (+) value indicates the element is in tension, while a (-) value in compression.\n")
    print('\u2588'*100)

