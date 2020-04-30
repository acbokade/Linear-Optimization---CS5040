import numpy as np

def simplex_phase1(N, M, A, B, C) : 

    arr=[]
    for i in range(M):
        arr.append(i)
    data = [0]*N
    global extreme_pt
    combinationUtil(arr, data, 0,M - 1, 0, N, N, M, A, B, C);
    return extreme_pt


def isvalid(data, N, M, A, B, C):
    global extreme_pt
    mat = []
    tempb = []
    for i in range (N):
        mat.append(A[data[i]])
        tempb.append(B[data[i]])
    mat = np.array(mat, dtype='float')
    tempb = np.array(tempb, dtype='float')

    print(mat)
    if abs(np.linalg.det(mat)-0)<1e-5:
        return (False)

    matinverse= np.linalg.inv(mat)

    sol = np.dot(matinverse,tempb)
    
    tempmul = np.dot(A,sol)
    print(tempmul)
    for i in range(M):
        if (B[i]-tempmul[i] < 0):
          return False
    extreme_pt = sol    
    return (True)
    

def combinationUtil(arr, data, start,  
                    end, index, r, N, M, A, B, C): 
                            
    # Current combination is ready  
    # to be printed, print it 
    global extreme_pt
    if len(extreme_pt)!=0:
        return
    
    if index == r: 
        print(data)
        print("end\n")
        isvalid(data, N, M, A, B, C)
       
        print()
        return 
                
                
            
    
  
    # replace index with all 
    # possible elements. The 
    # condition "end-i+1 >=  
    # r-index" makes sure that  
    # including one element at 
    # index will make a combination  
    # with remaining elements at  
    # remaining positions 
    i = start;  
    while(i <= end and end - i + 1 >= r - index): 
        data[index] = arr[i]; 
        combinationUtil(arr, data, i + 1,  
                        end, index + 1, r, N, M, A, B, C); 
        i += 1; 



def check_tight(N,M,A,B,X):
    tight=[]
    for i in range(M):
        temp = 0
        for j in range(N):
            temp = temp +A[i][j]*X[j]
        if temp ==B[i]:
            tight.append(True)
        else:
            tight.append(False)
            
    return tight       

def simplex_phase2(X,tight_constraints, N, M, A, B, C):
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
            tight_constraints[untight_index_map[index_replace]] = True # untight becomes tight
            tight_constraints[tight_index_map[index]] = False # tight becomes untight


    if flag==1 :
        # print("Optimal Cost is" + optimal_cost + " at :" + X)
        return X

def Simplex(N, M, A, B, C):

    # find initial feasible solution
    # X is initial point and tight_constraints is boolean vector indicating which constraints are tight 
    trivial=1
    for i in range(M):
        if B[i] < 0:
            trivial=0
    if trivial == 1:
        X=[0]*N
    else:
        X= simplex_phase1(N, M, A, B, C)
        if (len(X)==0):
            print("Original problem is infeasible")
            return
    
    
    tight_constraints = check_tight(N,M,A,B,X)
    
    optimal_pt =  simplex_phase2(X,tight_constraints, N, M, A, B, C)
    print (optimal_pt,C)
#     print("Optimal Cost is" + np.dot(C,optimal_pt) + " at :" + optimal_pt)
    





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
    print(A)
    print(B)
    print(C)
    # N=3
    # M=7
    # A = [[-1 ,0,0],[0,-1,0],[0,0,-1],[1,0,0],[0,1,0],[0,0,1],[1,1,1]]
    # B = [0,0,0,2,2,2,-5]
    # C= [1,1,-1]
    

 
    Simplex(N, M, A,B,C)