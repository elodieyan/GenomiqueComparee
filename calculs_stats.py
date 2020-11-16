#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 16:03:19 2020

@author: Chanco
"""

#### Modules

import pandas as pd
import pickle
import os
from math import *
from statistics import mean,stdev
import sys

#### functions

def basics_stats(workpath):
    
    rep = os.scandir(workpath)
    sizes = []
    
    for file in rep:
        
        df = pd.read_csv(file,sep='\t', header = 0)
        size = df.shape[0]
        sizes.append(size)
    
    variance = stdev(sizes)
    avg = mean(sizes)
    maxi = max(sizes)
    mini = min(sizes)    
    
    print('maximus size =',maxi )
    print('minimum size =',mini )
    print('mean size = ', avg)
    print('standart deviation = ',variance)

def comparisons(directory_1,directory_2):
    
    percents_bbh = []
    
    for file1 in os.scandir(directory_1):
        
        df1 = pd.read_csv(file1,sep = '\t', header = 0)
        name = file1.name.split('.csv')[0]
        file2 = directory_2 + name + '_BBH.csv'
        df2 = pd.read_csv(file2, sep = '\t', header = 0)
        percents_bbh.append((df2.shape[0]/df1.shape[0]) * 100)
    
    variance = stdev(percents_bbh)
    avg = mean(percents_bbh)
    maxi = max(percents_bbh)
    mini = min(percents_bbh)
    
    print('maximus % =',maxi )
    print('minimum % =',mini )
    print('mean % = ', avg)
    print('% standart deviation = ',variance)
        
        
    
    
    
####  Main

directory_paths = [str(sys.argv[1]),str(sys.argv[2])]
print(directory_paths)
print('analyzed repertory : ',directory_paths[0])
basics_stats(directory_paths[0])

print('\nanalyzed repertory : ',directory_paths[1])

basics_stats(directory_paths[1])

print('\nrepertory comparisons :')

comparisons(directory_paths[0],directory_paths[1])
