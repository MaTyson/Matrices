##################################################################
# Projeto MS512 - Análise Numérica                               #
# Prof Maicon R. Correa                                          #
# Fatorações QR via Gram-Schmidt e transformações de Householder #
#                                                                #
# Bruno S. T. Pagotto 154891                                     #
# João Pedro N. S. 176146                                        #
# Lucas G. Machado 172776                                        #
# Matheus L. Bernardi 184331                                     #
##################################################################

import numpy
import random
import sys
import time
import numpy.linalg as la
from scipy.linalg import invhilbert as invH


def hilbert(n):
    # função que gera uma matriz de Hilbert de dimensão nxn
    # input: dimensão n
    # output: matriz de Hilbert H
    
    H = numpy.zeros([n,n])
    for i in range(n):
        for j in range(n):
            H[i,j] += 1/(i + j + 1)

    return H 

def scalar(v, u):
    # função que calcula o produto escalar entre dois vetores
    # inputs: vetores v e u de dimensão n
    # output: produto escalar s
    
    n = v.size
    s = 0
    for i in range(n):
        s += u[i]*v[i]
    
    return s

def mprod(A, B):
    # função que calcula o produto entre duas matrizes
    # inputs: matrizes A e B
    # output: matriz C = AB
    
    m = A.shape[0]
    p = B.shape[1]
    C = numpy.zeros([m,p])
    
    for i in range(m):
        for j in range(p):
            C[i,j] += scalar(A[i,:], B[:,j])
    
    return C 

def vecmat(v,A):
    # função que calcula produto entre vetor e matriz
    # inputs: vetor v e matriz A
    # output: vetor x = v^tA
    
    m,n = A.shape
    x = numpy.zeros([n,])
    for i in range(n):
        x[i] += scalar(v,A[:,i])
    
    return x

def matvec(A,v):
    # função que calcula produto entre matriz e vetor
    # inputs: matriz A e vetor v
    # output: vetor x = Av
    
    m,n = A.shape
    x = numpy.zeros([m,])
    for i in range(m):
        x[i] += scalar(A[i,:],v)
    
    return x

def tensor(v,u):
    # função que calcula o produto externo entre dois vetores
    # inputs: vetores v e u
    # output: tensor T
    
    n = v.size
    T = numpy.zeros([n,n])
    
    for i in range(n):
        T[i,:] += v[i]*u
        
    return T

def sign(x):
    # função sinal
    # input: escalar x
    # output: sign(x)
    
    if(x >= 0):
        return 1
    else:
        return -1

def retrosub(U,c):
    # função que resolve um sistema linear via retrosubstituição
    # inputs: matriz triangular superior U e vetor c
    # output: vetor x solução do sistema Ux = c
    
    n = U.shape[1]
    x = numpy.zeros([n,])
    x[n-1] = c[n-1]/U[n-1,n-1]
    
    for i in range(n-2,-1,-1):
        x[i] = c[i]
        for j in range(i + 1, n ):
            x[i] = x[i] - U[i,j] * x[j]
        x[i] = x[i]/U[i,i]
    
    return x

def sismat(U,B):
    # função que resolve um sistema matricial UX = B via retrosubstituição
    # inputs: matrizes U e B, com U triangular superior
    # output: matriz X solução do sistema
    
    m, n = B.shape
    X = numpy.zeros([m,n])
    for i in range(n):
        X[:,i] += retrosub(U,B[:,i])
    
    return X

def sub(L,c):
    # função que resolve um sistema linear via substituição direta
    # inputs: matriz L triangular inferior e vetor c
    # output: vetor x solução do sistema Lx = c
    
    n = L.shape[1] -1
    x = numpy.zeros([n + 1,])
    x[0] = c[0]/L[0,0]
    
    for i in range(1, n+1):
        x[i] = c[i]
        for j in range(i):
            x[i] -= L[i,j] * x[j]
        x[i] /= L[i,i]
    
    return x

def gram(A):
    # função da fatoração QR via Gram-Schmidt
    # input: matriz A mxn
    # outputs: matrizes Q mxn com colunas ortonormais e R nxn triangular superior
    
    m, n = A.shape
    Q = numpy.zeros([m,n])
    R = numpy.zeros([m,m])
    
    for k in range(m):
        Q[:,k] = A[:,k]
        for i in range(k):
            Q[:,k] -= Q[:,i]*scalar(Q[:,i],Q[:,k])
        r = scalar(Q[:,k],Q[:,k])
        if(r):
            Q[:,k] = Q[:,k]/numpy.sqrt(r)
        else:
            return 0
    
    R = mprod(Q.T,A)
    
    return Q, R    

def householder(A):
    # função da fatoração QR via Householder
    # input: matriz A mxn
    # outputs: matrizes Q ortogonal mxm e R mxn "triangular superior"
    
    m, n = A.shape  
    Q = numpy.eye(m)
    
    for i in range(n-(m == n)):
        H = numpy.eye(m)
        x = A[i:,i]
        beta = max(numpy.abs(x))
        if( beta == 0):
            gamma = 0
        else:
            x = x/beta
            tau = sign(x[0]) * numpy.sqrt( scalar(x,x))
            x[0] += tau
            gamma = x[0]/tau
            x[1:] /= x[0]
            x[0] = 1
            tau = tau * beta
            
        H[i:,i:] = H[i:,i:] - tensor(gamma * x ,vecmat(x.T,H[i:,i:]))
        Q = mprod(Q,H)
        A = mprod(H,A)
        
    return  Q, A

def norma_m1(A):
    # função que calcula norma 1 de matriz
    # input: matriz A
    # output: norma 1 de A
    
    m,n = A.shape
    x = numpy.zeros([n,])
    for i in range(n):
        x[i] += (numpy.abs(A[:,i])).sum()
    
    return numpy.max(x)

def norma_f(A):
    # função que calcula norma de Frobenius de matriz
    # input: matriz A
    # output: norma_f(A)
    
    return numpy.sqrt(((numpy.abs(A))**2).sum())

def norma_v1(v):
    # função que calcula norma 1 de vetor
    # input: vetor v
    # output: norma 1 de v
    
    return (numpy.abs(v).sum())

def cond(A,InvA):
    # função que calcula o número de condicionamento de A
    # inputs: matrizes A e sua inversa InvA
    # output: cond(A)
    
    return (norma_m1(A)*norma_m1(InvA))

def cond_e1(A):
    # função que calcula uma estimativa para o condicionamento de A
    # input: matriz A
    # output: estimativa cond(A) >= norma(a_i)/norma(a_j), com a_i e a_j colunas de A
    
    m,n = A.shape
    x = numpy.arange(n)
    x = random.sample(list(x),2)
    i = x[0]
    j = x[1]
    
    return norma_v1(A[:,i])/norma_v1(A[:,j])

def cond_e2(Q,R):
    # função que calcula uma estimativa para o condicionamento de A
    # inputs: matrizes Q e R da fatoração A = QR
    # output: estimativa cond(A) >= norma(A)*norma(InvA*w)/norma(w), para um vetor w qualquer
    
    A = mprod(Q,R)
    m,n = A.shape
    r = numpy.max(numpy.abs(A))
    w = r*numpy.random.rand(n,1)
    
    return norma_m1(A)*norma_v1(retrosub(R,matvec(Q.T,w)))/norma_v1(w)

if '__name__' == '__main__':
    t = time.time()
    
    # matriz de Hilbert H nxn
    n = int(sys.argv[1])
    H = hilbert(n)

    # fatorações QR
    Q1, R1 = gram(H)
    Q2, R2 = householder(H)

    # cálculo da inversa de H
    X1 = sismat(R1,Q1.T)
    X2 = sismat(R2,Q2.T)
    X3 = invH(n)

    # norma de HX, para X é a inversa de H
    n1 = norma_m1(mprod(H,X1))
    n2 = norma_m1(mprod(H,X2))
    n3 = norma_m1(mprod(H,X3))

    # condicionamento de H
    k1 = cond(H,X1)
    k2 = cond(H,X2)
    k3 = la.cond(H)

    # estimativas para o condicionamento de H
    e1 = cond_e1(H)
    e2 = cond_e2(Q1,R1)
    e3 = cond_e2(Q2,R2)

    t = time.time()-t
    print(numpy.round(t,4),'seconds')

    print(n1, n2, n3)
    print(k1, k2, k3)
    print(e1, e2, e3)