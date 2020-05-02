        #####################################################################################
        #              CS5040 - Linear Optimization                                         #
        #            Simplex Algorithm ( Degenerate )                                       #    
        #                                                                                   #    
        #              Ajinkya Bokade - CS17BTECH11001                                      #
        #             Shivashish Suman - CS17BTECH11037                                     #
        #            Rushikesh Tammewar - CS17BTECH11041                                    #
        #                                                                                   #
        #####################################################################################                 


import numpy as np
# boolean global variable for degenerate case
isDegenerate = False

# boolean global variable for unbounded solution
isUnbounded = False


def simplex_phase1(N, M, A, B, C):
    arr = []                                                                #  Returns First Extreme Point                                                 
    for i in range(M):
        arr.append(i)
    data = [0]*N
    global extreme_pt
    combinationUtil(arr, data, 0, M - 1, 0, N, N, M, A, B, C)               #  recursively finds out which N constraints are 
    return extreme_pt                                                       #  independent and feasible.



def isvalid(data, N, M, A, B, C):                 # checks  weather given N constraints are indepndent and Feasible.
    global extreme_pt
    mat = []
    tempb = []
    for i in range(N):
        mat.append(A[data[i]])
        tempb.append(B[data[i]])
    mat = np.array(mat, dtype='float')
    tempb = np.array(tempb, dtype='float')

    if abs(np.linalg.det(mat)-0) < 1e-5:
        return (False)

    matinverse = np.linalg.inv(mat)
    sol = np.dot(matinverse, tempb)
    tempmul = np.dot(A, sol)

    for i in range(M):
        if (B[i]-tempmul[i] < 0):
            return False
    extreme_pt = sol                        # intersection of N independent Feasible constraint
    return (True)                           # is assigned as extreme point. 


def combinationUtil(arr, data, start,
                    end, index, r, N, M, A, B, C):     # All combinations N out of M constraint are verified recursively.
    global extreme_pt
    if len(extreme_pt) != 0:
        return

    if index == r:
        isvalid(data, N, M, A, B, C)
        return

    i = start
    while(i <= end and end - i + 1 >= r - index):
        data[index] = arr[i]
        combinationUtil(arr, data, i + 1,
                        end, index + 1, r, N, M, A, B, C)
        i += 1


def check_tight(N, M, A, B, X):          # returns boolean list of tight and untight constraint
    tight = []
    for i in range(M):
        temp = 0
        for j in range(N):
            temp = temp + A[i][j]*X[j]
        if temp == B[i]:
            tight.append(True)
        else:
            tight.append(False)
    return tight


def simplex_phase2(X, tight_constraints, N, M, A, B, C):
    optimal_cost = np.dot(C, X)
    flag = -1

    while flag == -1:
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
            else:
                a_untight.append(A[i])
                b_untight.append(B[i])
                untight_index_map[untight_ind] = i
                untight_ind += 1

        a_tight = np.array(a_tight, dtype='float')
        a_untight = np.array(a_untight, dtype='float')
        b_tight = np.array(b_tight, dtype='float')
        b_untight = np.array(b_untight, dtype='float')
        global extreme_pt
        global isDegenerate
        global isUnbounded
        # degenerate case
        if len(a_tight) > N:
            extreme_pt = []
            isDegenerate = True
            epsilon = 1e-5
            # perturbing b vector, adding infinitesimally small epsilon to b vector
            for i in range(len(B)):
                B[i] = B[i] + epsilon ** (i+1)
                # solving new problem using simplex since its non degenerate
            X = Simplex(N, M, A, B, C)
            # flag = 1
            break
        alpha = np.matmul(np.linalg.inv(np.transpose(a_tight)), C)
        t = 1e18
        index_replace = -1
        for i in range(N):
            if alpha[i] < 0:
                index = i
                break

        if index == -1:
            flag = 1
            break
        else:
            a_inverse = np.linalg.inv(a_tight)
            # direction which increases the cost
            dir_vect = -a_inverse[:, index]
            # finding the untight row to replace with tight row
            for i in range(b_untight.size):
                val = np.dot(a_untight[i], dir_vect)
                # if cost is increasing in direction vector and no constraint is becoming tight,
                # it means it is unbounded solution
                if val <= 0:
                    isUnbounded = True
                    break
                temp = (b_untight[i] - np.dot(a_untight[i], X))/val
                if t > temp:
                    t = temp
                    index_replace = i
                        
            if isUnbounded:
                print("Unbounded Solution")
                break
            # degenerate case when more than 1 row is becoming tight
            # that is more than n tight rows
            increase_in_tight_rows = 0
            for i in range(b_untight.size):
                val = np.dot(a_untight[i], dir_vect)
                temp = (b_untight[i] - np.dot(a_untight[i], X))/val
                # if t == temp:
                #     increase_in_tight_rows += 1
                if abs(t-temp) <= 1e-05:
                    increase_in_tight_rows += 1

            # degenerate case, perturbing the b vector
            if increase_in_tight_rows > 1:
                extreme_pt = []
                isDegenerate = True
                epsilon = 1e-5
                # perturbing b vector, adding infinitesimally small epsilon to b vector
                for i in range(len(B)):
                    B[i] = B[i] + epsilon ** (i+1)
                # solving new problem using simplex since its non degenrate
                X = Simplex(N, M, A, B, C)
                break

            X = X+dir_vect*t
            # untight becomes tight
            tight_constraints[untight_index_map[index_replace]] = True
            tight_constraints[tight_index_map[index]
                              ] = False  # tight becomes untight

    if flag == 1:
        if(isDegenerate):
            print("Degenerate case")
        print("Optimal Cost is " + str(np.dot(C, X)) +
              " at :" + str(X))
        return X


def Simplex(N, M, A, B, C):

    # find initial feasible solution
    # X is initial point and tight_constraints is boolean vector indicating which constraints are tight
    X = simplex_phase1(N, M, A, B, C)
    # if no n linearly independent tight constraints are found, then solution is infeasible
    if (len(X) == 0):
        print("Original problem is infeasible")
        return

    tight_constraints = check_tight(N, M, A, B, X)

    optimal_pt = simplex_phase2(X, tight_constraints, N, M, A, B, C)


if __name__ == "__main__":

    N = int(input("Enter the number of variables (n):"))
    M = int(input("Enter the number of constraints (m):"))

    print("Enter the matrix A row-wise:")
    A = [[j for j in raw_input().strip().split(" ")] for i in range(M)]

    print("Enter the elements of B separated by space:")
    B = [j for j in raw_input().strip().split(" ")]

    print("Enter the elements of C separated by space:")
    C = [j for j in raw_input().strip().split(" ")]

    A = np.array(A, dtype='float')
    B = np.array(B, dtype='float')
    C = np.array(C, dtype='float')
    extreme_pt = []
    Simplex(N, M, A, B, C)
