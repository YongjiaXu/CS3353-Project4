#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 17:35:11 2020

@author: yongjiaxu
"""
import numpy as np
import os.path
import random 
counter = 0

def generateMatrix(size):
    
    norm = 10
    c = 0
    while(norm >=1):
        v = []
        for i in range(size):
            temp = []
            for j in range(size):
                x = random.randint(1,10)
                temp.append(x)
            v.append(temp)
                        
        m = np.max(v)
        
        #print(np.array(v))
        for i in range(size):
            v[i][i] *= m
        
        #check if the matrix converges
        l = np.tril(v)
        u = np.triu(v)
        d = np.diag(v)
        
        # scalar for u in order for it to converge faster
        u = u * 0.01
        
        try:
            -np.dot(np.linalg.inv(d+l), u)
        except np.linalg.LinAlgError:
            continue
        
        b = -np.dot(np.linalg.inv(d+l), u)
        norm = np.linalg.norm(b)
        print(norm)
        c += 1
        
        v = l+d+u

    return v
        

def writeFile(size, counter, save_path = '/Users/yongjiaxu/Desktop/CS3353/Program4/temp'):
    
    completeName = os.path.join(save_path, "{}_{}.txt".format(size, counter))         
    
    file1 = open(completeName, "w")
    
    file1.write(str(size))
    file1.write('\n')

    v = generateMatrix(size)
    for i in range(size):
        for j in range(size):
            file1.write(str(v[i][j])+' ')
        x = random.random()
        file1.write(str(x))
        file1.write('\n')

    file1.close()
    
    
        
    
    
    
    
    
    
    
    
    
    
    
    