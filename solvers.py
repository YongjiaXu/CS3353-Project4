#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 01:57:15 2020

@author: yongjiaxu
"""
import numpy as np

def gaussElim(A):
    n = len(A)
    for i in range(0, n):
        # search for maximum in this column
        maxEl = abs(A[i,i])
        maxRow = i
        for k in range(i + 1, n):
            if abs(A[k,i]) > maxEl:
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
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = A[i,n] / A[i,i]
        for k in range(i - 1, -1, -1):
            A[k,n] -= A[k,i] * x[i]
            
    return np.transpose(x)


def gaussSeidel(a, x, b, tol = 1e-8):
    x1 = x.copy()
    for i in range(0, 100000):             
        x = inner(a, x, b) 
        err = max(np.absolute(x-x1))
        if(err < tol):
            break;
        x1 = x.copy()
        
    return x;


def inner(a, x ,b): 
    n = len(a)                    
    for j in range(0, n):         
        d = b[j]            
        for i in range(0, n):      
            if(j != i):
                # here x is still holding values from the previous solution
                d -= a[j,i] * x[i] 
        # updating the values of our solution         
        x[j] = d / a[j,j]
    return x