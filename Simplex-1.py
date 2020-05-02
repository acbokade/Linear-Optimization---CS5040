        #####################################################################################
        #              CS5040 - Linear Optimization                                         #
        #            Simplex Algorithm ( Non-degenerate )                                   #    
        #                                                                                   #    
        #              Ajinkya Bokade - CS17BTECH11001                                      #
        #             Shivashish Suman - CS17BTECH11037                                     #
        #            Rushikesh Tammewar - CS17BTECH11041                                    #
        #                                                                                   #
        #####################################################################################                                                                                   #    



import numpy as np

# boolean global variable for unbounded solution
isUnbounded = False     

def simplex_phase1(N, M, A, B, C) :                           #  Returns First Extreme Point               
    arr=[]
    for i in range(M):
        arr.append(i)
    data = [0]*N
    global extreme_pt
    combinationUtil(arr, data, 0,M - 1, 0, N, N, M, A, B, C)   #  recursively finds out which N constraints are 
    return extreme_pt                                          #  independent and feasible.


def isvalid(data, N, M, A, B, C):        # checks  weather given N constraints are indepndent and Feasible.
    global extreme_pt
    mat = []
    tempb = []
    for i in range (N):
        mat.append(A[data[i]])
        tempb.append(B[data[i]])
    mat = np.array(mat, dtype='float')
    tempb = np.array(tempb, dtype='float')

    if abs(np.linalg.det(mat)-0)<1e-5:
        return (False)

    matinverse= np.linalg.inv(mat)
    sol = np.dot(matinverse,tempb)
    
    tempmul = np.dot(A,sol)
    for i in range(M):
        if (B[i]-tempmul[i] < 0):
            return False
    extreme_pt = sol                  # intersection of N independent Feasible constraint    
    return (True)                     # is assigned as extreme point. 
    

def combinationUtil(arr, data, start,                  # All combinations N out of M constraint are verified recursively.
                    end, index, r, N, M, A, B, C): 
                            
    global extreme_pt
    if len(extreme_pt)!=0:
        return
    
    if index == r: 
        isvalid(data, N, M, A, B, C)
        return 
    
    i = start
    while(i <= end and end - i + 1 >= r - index): 
        data[index] = arr[i]
        combinationUtil(arr, data, i + 1,  end, index + 1, r, N, M, A, B, C)
        i += 1



def check_tight(N,M,A,B,X):    # returns boolean list of tight and untight constraint
    tight=[]
    for i in range(M):
        temp = 0
        for j in range(N):
            temp = temp +A[i][j]*X[j]
        if abs(temp - B[i]) < 1e-5:
            tight.append(True)
        else:
            tight.append(False)
            
    return tight       

def simplex_phase2(X,tight_constraints, N, M, A, B, C):
    flag = -1

    while flag == -1 :
        index = -1
        a_tight = []        # collection of tight rows
        b_tight = []        # collection of values of b corresponding to tight rows
        a_untight = []      # collection of untight rows
        b_untight = []      # collection of values of b corresponding to untight rows
        tight_index_map = {}   # map to store which rows are tight
        untight_index_map = {}   # map to store which rows are untight
        tight_ind = 0
        untight_ind = 0
        for i in range(M):
            if tight_constraints[i]:
                a_tight.append(A[i])
                b_tight.append(B[i])
                tight_index_map[tight_ind] = i
                tight_ind += 1
            else :
                a_untight.append(A[i])
                b_untight.append(B[i])
                untight_index_map[untight_ind] = i
                untight_ind +=1
        
        a_tight = np.array(a_tight, dtype='float')
        a_untight = np.array(a_untight, dtype='float')
        b_tight = np.array(b_tight, dtype='float')
        b_untight = np.array(b_untight, dtype='float')
        global isUnbounded
        
        alpha = np.matmul(np.linalg.inv(np.transpose(a_tight)),C)
        t = 1e18
        index_replace=-1
        for i in range(N):
            if alpha[i]<0:
                index=i
                break
        
        if index == -1 :
            flag=1
            break
        else :
            a_inverse = np.linalg.inv(a_tight)
            dir_vect = -a_inverse[:,index]          # direction which increases the cost

            #finding the untight row to replace with tight row 
            for i in range(b_untight.size):
                val = np.dot(a_untight[i],dir_vect)
                if val <= 0:
                    continue
                temp = (b_untight[i] - np.dot(a_untight[i],X))/val 
                if t>temp:
                    t=temp
                    index_replace=i
            
            # if no untight row is becoming tight and cost is increasing
            # in the direction of direction vector, then it has unbounded solution
            if t == 1e18:
                isUnbounded = True
                print("Unbounded Solution")
                break    
                   
            X=X+dir_vect*t
            tight_constraints[untight_index_map[index_replace]] = True # untight becomes tight
            tight_constraints[tight_index_map[index]] = False # tight becomes untight


    if flag==1 :
        return X

def Simplex(N, M, A, B, C):

    # find initial feasible solution
    # X is initial Extreme point and tight_constraints is boolean vector indicating which constraints are tight 
    X= simplex_phase1(N, M, A, B, C) 
    # if no n linearly independent tight constraints are found, then solution is infeasible
    if (len(X)==0):
        print("Original problem is infeasible")
        return
    tight_constraints = check_tight(N,M,A,B,X)
    
    #Simplex phase 2 returns optimal point 
    optimal_pt =  simplex_phase2(X,tight_constraints, N, M, A, B, C)

    print("Optimal Cost is " + str(np.dot(C,optimal_pt)) + " at :" )
    print(optimal_pt)
    





if __name__ == "__main__":

    N = int(input("Enter the number of variables (n):")) 
    M = int(input("Enter the number of constraints (m):")) 

    print("Enter the matrix A row-wise:") 
    A = [[j for j in input().strip().split(" ")] for i in range(M)] 

    print("Enter the elements of B separated by space:")
    B = [j for j in input().strip().split(" ")]

    print("Enter the elements of C separated by space:")
    C = [j for j in input().strip().split(" ")]

    A = np.array(A, dtype='float')
    B = np.array(B, dtype='float')
    C = np.array(C, dtype='float')

    extreme_pt = []
    
    Simplex(N, M, A,B,C)
