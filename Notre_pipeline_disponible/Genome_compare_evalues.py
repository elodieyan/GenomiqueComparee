#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 11:33:12 2020

@author: Chanco
"""

### Imports

import pandas as pd
import pickle 
import progressbar
from io import StringIO
import sys
import os

### 

def Best_hit(path_to_file, ev_treshold, ):
    
    file = open(path_to_file,'r')
    txt = file.read()
    querys = txt.split('# BLASTP 2.2.31+')[1:-1]
    string =''
    string += 'query id\tsubject id\t% identity\talignment length\tmismatches\tgap opens\tgaps\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\tquery length\tsubject length' + '\n'
    for q in querys:
        lines = q.split('\n')[1:]
        if lines[2] == '# 0 hits found':
            continue
        else:
            best_hit = lines[4]
            string += best_hit + '\n'
    
    dt = StringIO(string)

    df = pd.read_csv(dt, header = 0,sep="\t")
    df['evalue'] = df['evalue'].astype(float)
    ndf = df[['query id','subject id','evalue']]
    ### ICI on doit faire la selection cds vs cds , igor vs igor bref
    
    
    return ndf