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
    files = [f for f in glob.glob(path + "**/*3_0.txt", recursive=True)]
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
#    general_test()
    
    for i in [0,1,2,3,4]:
        start = timer()
        sol,it = ode_2nd.solve(i, 5)
        end = timer()
        print(i, sol,end - start)
    
    #generate converging matrix for gauss seidel
#    for i in [3,4,5]:
#        print(i)
#        matGenerator.writeFile(i, save_path = '/Users/yongjiaxu/Desktop')
    
        

    
        
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
