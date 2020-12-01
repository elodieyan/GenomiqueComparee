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

def Best_hit(path_to_file):
    
    file = open(path_to_file,'r')
    txt = file.read()
    querys = txt.split('# BLASTP 2.6.0+')[1:-1]
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
    
    nature_query = [len(i.split('_')) for i in list(ndf['query id'])]
    nature_subject = [[len(i.split('_')) for i in list(ndf['subject id'])]]
    
    ndf['nature_query'] = nature_query
    ndf['nature_subject'] = nature_subject
    
    ndf_cds_vs_cds = ndf[(ndf['nature_query'] == 3) & (ndf['nature_subject'] == 3)]
    ndf_igorf_vs_igorf = ndf[(ndf['nature_query'] == 4) & (ndf['nature_subject'] == 4)]             
    
    return ndf