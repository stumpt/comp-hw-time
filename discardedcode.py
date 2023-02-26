j = 0
while j < numNodes*2:
    print(str(j))
    if j in nodesWithLoads:
        j+=2
    #first case: displacement BC = 0. In this case, set combined force values to 0, and alter 
    #globalK to reflect the changes.
    elif uvList[j] == 0:
        print("u = 0 at row: " + str(j))
        forceList[j] = 0 
        globalK[j] = 0
        globalK[j,j] = 1
        j+=1
    elif uvList[j] == None:
        print("none detected at row " + str(j))
        j+=1
    else: #when we have a known, non-zero displacement boundary condition
        #alter global stiffness to all zeros and one 1, and set force equal to displacement BC
        print("row " + str(j) + " altered")
        print(uvList[j])
        forceList[j] = uvList[j]
        for c in range(numNodes):
            if (globalK[j,c] > 0) == True:
                globalK[j,c] = 1
            else:
                globalK[j,c] = 0
        j+=1

                
print(globalK)
print(forceList)