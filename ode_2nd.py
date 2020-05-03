#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  2 19:00:01 2020

@author: yongjiaxu
"""

import numpy as np
from numpy import exp,zeros,arange,absolute,max
from math import exp,log,sqrt,cos,pi,sin,pi
import matplotlib.pyplot as plt
from timeit import default_timer as timer
#from ODEs import *   #import the self-coded ODE solvers
from scipy.integrate import  quad, odeint

# solution function
def func(x):
    y= (2*exp(1))*x*(exp(-x)) - exp(x)
    return y

# righthand side function
def rfunc(x):
    y= -4*(exp(x))
    return y

#y = [y y']
def f(y,x):
    return np.array([y[1], -4*np.exp(x) - 2*y[1] - y[0]])

def solve(method, n):
#    print(n)'
    xa=0; xb=2;
    h = (xb-xa)/(n+1);  # mesh size
    alpha = func(xa); beta = func(xb)
    x = np.linspace(xa, xb, n)
    if (method == 5):
        s = odeint(f, alpha, x)
        sol = s[:,0]
        return sol, -1
    A=1; B=2; C=1

    # exact value at a fine mesh
#    d=0.0025
#    xe = arange(xa,xb+d,d)
#    ye = xe.copy()
#    for i in range(len(xe)):
#        ye[i] = func(xe[i])

    # matrix entry on tri-dialgonals
    coA= A/(h**2) - B/(2*h)
    coB= C-2*A/(h**2)
    coC= A/(h**2) + B/(2*h)

    xh = zeros(n) #x-values
#    yh = zeros(n) #true y-values at grids
    sol = zeros(n) #computed y-values at grids
    r=zeros(n)    #right handside of the equations

    for i in range(n):
        xh[i] = (i+1) * h + xa
#        yh[i] = func(xh[i])
        r[i] = rfunc(xh[i])
    r[0] = r[0] - alpha*coA
    r[n-1] = r[n-1] - beta *coC

    if (method == 0):
        # Thomas's algorithm
        # assign values for a,b,c - tridiagonals
        a=zeros(n);b=zeros(n);c=zeros(n)
        a += coA; b += coB; c += coC
        sol = Thomas(a,b,c,r)
        iteration = -1
    elif(method == 1):
        # Jacobi iterative method
        sol, iteration = Jacobi(coA, coB, coC, r, n)
    elif(method == 2):
        # Gauss-Seidel iterative method
        sol, iteration = Gauss_Seidel(coA, coB, coC, r, n)
    elif(method == 3):
        # Successive over-relaxation
        wopt = 2/(1+sqrt(1-cos(pi*h)**2)) #optimized omega for tridiagonal system
        sol, iteration = SOR(coA, coB, coC, r, n, wopt)
    elif(method == 4):
        # Gaussian Elimination:
        M = np.zeros((n,n+1))
        for i in range(1,n-1):
            M[i,i] = coB
            M[i, i-1] = coA
            M[i, i+1] = coC
        M[0,0] = coB
        M[0, 1] = coC
        M[n-1, n-2] = coA
        M[n-1,n-1] = coB
        for i in range(n):
            M[i,n] = r[i]
        sol = gaussElim(M)
        iteration = -1
    return sol, iteration

def Jacobi(coA, coB, coC, r, n, tol = 1e-8):
    # sol1 and sol are two consequent solutions
    # if the difference between two is smaller than tolerance 
    # then assume the result has been found
    sol1 = np.zeros(n)
    sol = np.zeros(n)
    iteration = 0
    err = 1
    while (err>tol):
        iteration += 1
        #calculation part
        sol[0] = (r[0] - coC * sol1[1])/coB
        for i in range(1, n-1):
            sol[i] = (r[i] - coA * sol1[i-1] - coC * sol1[i+1])/coB
        sol[n-1] = (r[n-1] - coA * sol1[n-2])/coB
        
        err = max(absolute(sol1-sol))
        sol1 = sol.copy()
        
    return sol, iteration


def Gauss_Seidel(coA, coB, coC, r, n, tol = 1e-8):
    # sol1 and sol are two consequent solutions
    # if the difference between two is smaller than tolerance 
    # then assume the result has been found
    sol1 = np.zeros(n)
    sol = np.zeros(n)
    iteration = 0
    err = 1
    while (err>tol):
        iteration += 1
        #calculation part 
        #** notice here we use only sol because the only place we use the previous solution
        # is to use the next value and since that value in sol is not updated
        # we could just use the same array
        sol[0] = (r[0] - coC * sol[1])/coB
        for i in range(1, n-1):
            sol[i] = (r[i] - coA * sol[i-1] - coC * sol[i+1])/coB
        sol[n-1] = (r[n-1] - coA * sol[n-2])/coB
        
        err = max(absolute(sol1-sol))
        sol1 = sol.copy()
        
    return sol, iteration

def SOR(coA, coB, coC, r, n, wopt, tol = 1e-8):
    # sol1 and sol are two consequent solutions
    # if the difference between two is smaller than tolerance 
    # then assume the result has been found
    sol1 = np.zeros(n)
    sol = np.zeros(n)
    iteration = 0
    err = 1
    while (err>tol):
        iteration += 1
        #calculation part 
        #** notice here we use only sol because the only place we use the previous solution
        # is to use the next value and since that value in sol is not updated
        # we could just use the same array
        sol[0] = (r[0] - coC * sol[1])/coB
        for i in range(1,n-1):
            sol[i] = (coB*sol[i] + (wopt * (r[i] - coA * sol[i-1] - coB * sol[i] - coC * sol[i+1])))/coB
        sol[n-1] = (r[n-1] - coA * sol[n-2])/coB
        
        err = max(absolute(sol1-sol))
        sol1 = sol.copy()
        
    return sol, iteration 
    

#LU factorization based on Thomas's algorithm
def Thomas(a,b,c,r):
    n = len(r)
    sol = np.zeros(n)
    for k in range(1, n):
        q = a[k]/b[k-1]
        b[k] = b[k] - c[k-1]*q
        r[k] = r[k] - r[k-1]*q
    
    q = r[n-1]/b[n-1]
    sol[n-1] = q
    
    for k in range(n-2, -1, -1):
        q = (r[k]-c[k] * q)/ b[k]
        sol[k] = q
    return sol


# adapted to solve numpy array matrix
def gaussElim(A):
    n = len(A)
    for i in range(0, n):
        # search for maximum in this column
        maxEl = abs(A[i,i])
        maxRow = i
        for k in range(i + 1, n):
            if abs(A[k][i]) > maxEl:
                maxEl = abs(A[k,i])
                maxRow = k
                
        # check if there is any zero on the diagonal
        if A[i,maxRow] == 0:
            print('Singular')
            return None

        # swap maximum row with current row (column by column)
        for k in range(i, n + 1):
            tmp = A[maxRow,k]
            A[maxRow,k] = A[i,k]
            A[i,k] = tmp

        # make all rows below this one 0 in current column
        for k in range(i + 1, n):
            c = -A[k,i] / A[i,i]
            for j in range(i, n + 1):
                if i == j:
                    A[k,j] = 0
                else:
                    A[k,j] += c * A[i,j]

    # Solve equation Ax=b for an upper triangular matrix A
    x = [0 for i in range(n)]
    for i in range(n - 1, -1, -1):
        x[i] = A[i,n] / A[i,i]
        for k in range(i - 1, -1, -1):
            A[k,n] -= A[k,i] * x[i]
    return np.transpose(x)


if __name__ == "__main__":
#    solve_demo()
    print('\n')
    for i in [0,1,2,3,4]:
        start = timer()
        sol,it =solve(i, 5)
        end = timer()
        print(i, sol,end - start)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    