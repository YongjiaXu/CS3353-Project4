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
import matGenerator
import glob
import solvers
import ode_2nd

def general_test(path = '/Users/yongjiaxu/Desktop/CS3353/Program4/CS3353-Project4/matrices0'):
    s = []; GEtime = {}; GStime = {}; scipyTime = {}   
    files = [f for f in glob.glob(path + "**/*.txt", recursive=True)]
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
        solvers.gaussElim(A)
        end1 = timer()
        GEt = end1 - start1
        
        start2 = timer()
        solvers.gaussSeidel(a, x, b)
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
#    plt.savefig('GEvsGSvsScipy.png')
    plt.show()
    return



if __name__ == "__main__":
    print('-------------------------------------------------------------------------------')
    print('1. Run generate test (Gaussian Elimination vs Gauss-Seidel)')
    print('2. Run special numerical problem - second order ordinary differential equation')
    print('-------------------------------------------------------------------------------')
    option = input('Your option: ')
    
    if(option == '1'):
        # Part I: general matrix solving
        print('Do you have a particular path to your input matrices?')
        while(True):
            yn = input('Yes or No (Y/N):')
            if (yn == 'y') | (yn == 'Y') | (yn == 'yes') | (yn == 'Yes'):
                path = input('Please enter your directory: ')
                general_test(path)
                break;
            elif (yn == 'n') | (yn == 'N') | (yn == 'no') | (yn == 'No'):
                print('Using default path...')
                general_test()
                break;
        
    elif (option == '2'):
        # Part II: a specific numerical problem: 2nd order differential equation
        # Ay'' + By' + C = r(x); y(a) = alpha; y(b) = beta; 
        # Requirement: analytical solution: s1; right hand side: r(x); coefficients [A B C];
        #  endpoints xa, xb
        
        
        s1 = lambda x: (2*np.exp(1))*x*(np.exp(-x)) - np.exp(x)
        rf1 = lambda x: -4*(np.exp(x))
        xa1 = 0; xb1 = 2; coeff1 = [1,2,1]
        
        
        # ode_interface requires arguments: initial values of the ode: s, rf, xa, xb, coeff
        # it also has default arguments: n, m, T, I, E
        # n - data points; m - methods; T - display runtim; I - display # of iterations; E - error
        # default settings are n = [10,20,30]
        ode_2nd.ode_interface(s1, rf1, xa1, xb1, coeff1, n = [150,175,200,225,250], T = False, I = False)

    
    #generate converging matrix for gauss seidel
#    for i in [3,4,5]:
#        print(i)
#        matGenerator.writeFile(i, save_path = '/Users/yongjiaxu/Desktop')
    
        

    
        
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
