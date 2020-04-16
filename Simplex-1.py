import numpy as np

def find_initial(N, M, A, B, C) : 
    return A,B

def Simplex(N, M, A, B, C):

    # find initial feasible solution
    X, Basis = find_initial(A, B, C)

    Z = []

    for i in range(N):
        if Basis[i]:
            Z.append(A[:,i])
    Z = np.transpose(Z)
    flag = -1

    while flag == -1 :
        index = -1
        for i in range(N):
            if not Basis[i]:
                c_prime = c[i] - np.matmul(np.linalg.inv(Z),A[:,i])
                if c_prime 




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