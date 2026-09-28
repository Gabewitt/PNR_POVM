import numpy as np
import math
import scipy.special


M = 28 #pixel number
eta = 0.85 #detector efficiency
m = M #photon sent
precision = 4
np.set_printoptions(precision)
np.set_printoptions(suppress=True)
np.set_printoptions(linewidth=(precision+10)*M+5)
def P_zero(eta, n):
    return (1-eta)**n

def P_nn(eta, M, n):
    return (eta/M)**n * scipy.special.factorial(M)/scipy.special.factorial(M-n)

def povm_recursion_time(eta, M, m):
    P = np.zeros((M+1,m+1))

    for i in range(M+1):
        P[0][i] = P_zero(eta, i)

    for i in range(1, M+1):
        P[i][i] = P_nn(eta, M, i)

    for i in range(1,M+1): #photon detected
        for j in range(i+1,m+1): #photon sent
            P[i][j] = P[i][j-1] * ( (1-eta) + eta * i / M ) + P[i-1][j-1] * ( (M+1-i) * eta/M)
    return P


def probability(eta, MM, m, n):
    """MM is the number of pixel"""
    """n photon sent"""
    """m photon detected"""
    """eta : efficiency of the array"""
    p = 0
    for j in range(m+1):
        p += (-1)**j * ( scipy.special.factorial(m) / ( scipy.special.factorial(j) * scipy.special.factorial(m-j) ) ) * ( 1- eta + (m-j)*eta/MM)**n
    pp = scipy.special.factorial(MM) / ( scipy.special.factorial(m) * scipy.special.factorial (MM-m)) * p
    return pp


def povm_uniform_distibution(eta, M):
    P = np.zeros((M+1,M+1))
    for m in range(0,M+1): #photon I detect
        p = []
        for n in range(m,M+1): #photon I sent
            pp = probability(eta,M,m,n)
            P[m][n]= pp
    return P


# """Formula from Matioli et al 2009--> not working"""
# def povm_recursion_pixel(eta, M, m): 
#     P = np.zeros((M+1,m+1))

#     for i in range(M+1):
#         P[0][i] = P_zero(eta, i)

#     for i in range(1, M+1):
#         P[i][i] = P_nn(eta, M, i)

#     for i in range(1,M+1): #photon detected
#         for j in range(i+1,m+1): #photon sent
#             P[i][j] = P[i][j-1] * ( i/M + scipy.special.factorial(i)/scipy.special.factorial(M) * (1-eta/M)/M * (M-i) ) + P[i-1][j-1] * ( scipy.special.factorial(i-1)/scipy.special.factorial(M) * (eta/M)*(M-i+1))
#     return P


if __name__ == "__main__":
    eta = 0.90
    M = 28
    P = povm_uniform_distibution(eta, M)
    P_42 = povm_uniform_distibution(eta, 42)

    tes = []
    for i in range(10):
        tes.append(eta**i)
    # print(P)


    # for i in range(len(P[0])):
    #     c = P[:,i]
    #     print(sum(c))

    # name = f"POVM_{M}pixels_{eta}"
    # complete_name = name + ".txt"
    # file=open(complete_name,'w')
    # for i in range(len(P)):
    #     file.write(str(P[i]))
    #     file.write("\n")
    # file.close()
    for i in range(10):
        print(f"P_{i,i} is {P[i][i]*100}")
        print(f"P_{i,i} is {P_42[i][i]*100}")
        print(f"P_{i,i} is {tes[i]*100}")
        print("\n")