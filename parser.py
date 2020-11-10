#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 10:16:05 2020

@author: Chanco
"""

#### pour l'utiliser:

"""
Lancer dans le terminal:
    
python3 path/to/parser.py path/to/file/to/parse path/to/result/file

"""


#### imports

import sys
import pandas as pd
from io import StringIO

####

#file = open(str(sys.argv[1]),'r')
#bests_hit_file = open(str(sys.argv[2]),'w')
file = open('/Users/Chanco/Desktop/AgroParisTech/3A/AMI2B/Genomique_comparee/full_pipeline/Blast_output/Escherichia_coli_536-vs-Escherichia_coli_55989.txt','r')

txt = file.read()
querys = txt.split('# BLASTP 2.11.0+')[1:-1]
database_info = querys[0].split('\n')[1:][1]
string =''
string += database_info + '\n'
string += 'query id\tsubject id\t% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\tquery length\tsubject length\tgaps' + '\n'



for q in querys:
    lines = q.split('\n')[1:]
    if lines[2] == '# 0 hits found':
        #best_hit = lines[0][9:]+ '\t' +'no hits found'
        #querys_with_no_hits_found.append(lines[0][9:])
        continue
    else:
        best_hit = lines[4]
        string += best_hit + '\n'
    
### Critères de séléction

dt = StringIO(string)

df = pd.read_csv(dt, header = 1,sep="\t")
df['% couverture'] = (df['alignment length'] / max([df['query length'],df['subject length']]) *100
ndf = df[df['% identity'] > 70 and df['evalue'] < 0.01 and df['% couverture'] >  60 ]

print(df[df['% identity'] > 70])


"""
% identity : 70%
evalue : 0,01
à définir seuil sur alignment length


"""