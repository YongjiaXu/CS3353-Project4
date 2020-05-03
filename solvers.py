#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 16:35:45 2020

@author: yongjiaxu
"""

import matplotlib.pyplot as plt
import scipy.linalg as linalg 
import numpy as np
from timeit import default_timer as timer
import glob
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
            
    return x


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


def general_test():
    path = '/Users/yongjiaxu/Desktop/CS3353/Program4/matrices'
    GEtime = {}; GStime = {}; scipyTime = {}
    s = []
    files = [f for f in glob.glob(path + "**/*_0.txt", recursive=True)]
    for f in files:
        f = open(f)
        # extract matrix from the file
        size = int(f.readline())
        s.append(size)
        A = []; a = []; b = []
        for i in range(size):
            line = f.readline()
            lineList = line.split(" ")
            temp = []
            for k in lineList:
                temp.append(float(k))
            A.append(temp)
            temp2 = temp.copy()
            b.append(temp2.pop())
            a.append(temp2)
      
        #convert to ndarray for unification
        A = np.array(A)
        a = np.array(a)
        b = np.array(b)
        x = np.zeros(size)

        start1 = timer()
        gaussElim(A)
        end1 = timer()
        GEt = end1 - start1
        
        start2 = timer()
        gaussSeidel(a, x, b)
        end2 = timer()
        GSt = end2 - start2
        
        start3 = timer()
        c,_,_,_ = linalg.lstsq(a, b)
        end3 = timer()
        scipyT = end3 - start3 
        
        # since the sequence of files are random; dictionary keeps track of the size & timing
        GEtime[size] = GEt
        GStime[size] = GSt
        scipyTime[size] = scipyT
        
#        print ('Gauss Elimination: ', GEtime)
#        print ('Gauss-Seidel: ',GStime)
#        print ('Scipy: ',scipyTime)        
        
        f.close()
        
    s.sort()
    GEy = []
    GSy = []
    scipyy = []
    for i in s:
        GEy.append(GEtime[i])
        GSy.append(GStime[i])
        scipyy.append(scipyTime[i])

    # plot the data
    plt.plot(s, GEy, '-*', markersize = 5, label = 'GE')
    plt.plot(s, GSy, '-x', markersize = 5, label = 'GS')
    plt.plot(s, scipyy, '-o', markersize = 3, label = 'scipy')
    plt.legend()
    plt.xlabel('size of matrix (nxn)')
    plt.ylabel('runtime (s)')
    plt.savefig('GEvsGSvsScipy.png')
    plt.show()
    return



if __name__ == "__main__":
    general_test()
        

    
        
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
