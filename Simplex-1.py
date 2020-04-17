import numpy as np

def find_initial(N, M, A, B, C) : 
    # fill the code here and remove the return below
    return A,B 

def Simplex(N, M, A, B, C):

    # find initial feasible solution
    # X is initial point and tight_constraints is boolean vector indicating which constraints are tight
    X, tight_constaints = find_initial(N, M, A, B, C)
    optimal_cost = np.dot(C,X)
    flag = -1

    while flag == -1 :
        index = -1
        a_tight = []
        b_tight = []
        a_untight = []
        b_untight = []
        tight_index_map = {}
        untight_index_map = {}
        tight_ind = 0
        untight_ind = 0
        for i in range(M):
            if tight_constaints[i]:
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
            dir_vect = -a_inverse[:,index]
            for i in range(b_untight.size):
                val = np.dot(a_untight[i],dir_vect)
                if val<=0:
                    continue
                temp = (b_untight[i] - np.dot(a_untight[i],X))/val
                if t>temp:
                    t=temp
                    index_replace=i
            X=X+dir_vect*t
            tight_constaints[untight_index_map[index_replace]] = True # untight becomes tight
            tight_constaints[tight_index_map[index]] = False # tight becomes untight


    if flag==1 :
        print("Optimal Cost is" + optimal_cost + " at :" + X)




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
    print(A)
    print(B)
    print(C)
    Simplex(N, M, A,B,C)