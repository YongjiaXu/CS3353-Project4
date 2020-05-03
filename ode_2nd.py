#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  2 19:00:01 2020

@author: yongjiaxu
"""

import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer
#from ODEs import *   #import the self-coded ODE solvers
from scipy.integrate import odeint
import solvers

def solve(method, n, rfunc, xa, xb, A, B, C, alpha, beta):
#    print(n)'
    h = (xb-xa)/(n+1);  # mesh size
        
    # matrix entry on tri-dialgonals
    coA= A/(h**2) - B/(2*h)
    coB= C-2*A/(h**2)
    coC= A/(h**2) + B/(2*h)

    x = np.zeros(n) #x-values
    sol = np.zeros(n) #computed y-values at grids
    r = np.zeros(n)    #right handside of the equations

    for i in range(n):
        x[i] = (i+1) * h + xa
        r[i] = rfunc(x[i])
    r[0] = r[0] - alpha*coA
    r[n-1] = r[n-1] - beta *coC
    
    if (method == 0):
        # Thomas's algorithm
        # assign values for a,b,c - tridiagonals
        a = np.zeros(n); b = np.zeros(n); c = np.zeros(n)
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
        wopt = 2/(1+np.sqrt(1-np.cos(np.pi*h)**2)) #optimized omega for tridiagonal system
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
        sol = solvers.gaussElim(M)
        iteration = -1
        

    return x, sol, iteration


#LU factorization based on Thomas's algorithm
def Thomas(a,b,c,r):
    n = len(r)
    w = np.zeros(n,float)
    l = np.zeros(n,float)
    u = np.zeros(n,float)
    z = np.zeros(n,float)
    #Determine L,U factors
    u[0] = b[0]
    
    for k in range(1, n):
        l[k] = a[k]/u[k-1]
        u[k] = b[k] - l[k]*c[k-1]
        
    # Solve Lz = r
    z[0] = r[0]
    for k in range(1, n):
        z[k] = r[k] - l[k]*z[k-1]

    # Solve Uw = z.
    w[n-1] = z[n-1]/u[n-1]
    for k in range(n-2,-1,-1):
        w[k] = (z[k] - c[k]*w[k+1])/u[k]
        
        
    return w


def Jacobi(coA, coB, coC, r, n, tol = 1e-8):
    # sol1 and sol are two consequent solutions
    # if the difference between two is smaller than tolerance 
    # then assume the result has been found
    sol1 = np.zeros(n)
    sol = np.zeros(n)
    iteration = 0
    err = 1
    while (err > tol):
        iteration += 1
        #calculation part
        sol[0] = (r[0] - coC * sol1[1])/coB
        for i in range(1, n-1):
            sol[i] = (r[i] - coA * sol1[i-1] - coC * sol1[i+1])/coB
        sol[n-1] = (r[n-1] - coA * sol1[n-2])/coB
        
        err = max(np.absolute(sol1-sol))
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
        
        err = max(np.absolute(sol1-sol))
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
        
        err = max(np.absolute(sol1-sol))
        sol1 = sol.copy()
        
    return sol, iteration 
    

def ode_perf(f1, rf1, xa, xb, coeff, n, m, T = True, I = True):

    alpha = f1(xa); beta = f1(xb)
    
#    d=0.0025; 
#    xe = np.arange(xa,xb+d,d)
#    ye = xe.copy()
#    for i in range(len(xe)):
#        ye[i] = f1(xe[i])
    l = []
    if T == True:
        time = {}; 
        for i in m:
            time[i] = l.copy()
    if I == True:
        iteration = {};
        for i in m:
            iteration[i] = l.copy()

    for k in range(len(n)):
        for i in m:
            start = timer()
            x, sol,it = solve(i, n[k], rf1, xa, xb, coeff[0], coeff[1], coeff[2], alpha, beta)
            end = timer()
#            y = np.zeros(k) #true y-values at grids
#            y = f1(x)
#            err = max(np.abs(y-sol))
            if T == True: time[i].append(end-start)
            if I == True: iteration[i].append(it)
        
    name = {}
    name[0] = 'Thomas'
    name[1] = 'Jacobi'
    name[2] = 'Gauss-Seidel'
    name[3] = 'SOR'
    name[4] = 'Gauss-Elimination'
    
    linestyle = {}
    linestyle[0] = 'r-'
    linestyle[1] = 'b-'
    linestyle[2] = 'g-'
    linestyle[3] = 'y-'
    linestyle[4] = 'c-'
    if T == True:
        for i in m:
            plt.plot(n, time[i], linestyle[i], label = name[i])
        plt.legend()
        plt.title('Runtime performance')
        plt.xlabel('# of data points')
        plt.ylabel('Runtime(s)')
        plt.show()

    if I == True:
        for i in m:
            plt.plot(n, iteration[i], linestyle[i], label = name[i])  
        plt.title('# of iterations comparison')
        plt.xlabel('# of data points')
        plt.ylabel('# of iterations')
        plt.legend()
        plt.show()
    
    
def ode_interface():
    # f1 is the actual solution to the ode
    f1 = lambda x: (2*np.exp(1))*x*(np.exp(-x)) - np.exp(x)
    rf1 = lambda x: -4*(np.exp(x))
    n = [10,20,30]
    m = [0,3,4]
    coeff = [1,2,1]
    ode_perf(f1, rf1, 0, 2, coeff, n, m, T = True, I = False)

ode_interface()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    