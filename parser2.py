# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 10:08:13 2020

@author: Administrator
"""

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
file = open(r'C:\Users\Administrator\Desktop\AMI2B\Gen_Comp\blast_outputs\Escherichia_coli_536-vs-Escherichia_coli_55989.txt','r')

txt = file.read()
querys = txt.split('# BLASTP 2.2.31+')[1:-1]
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
df['% couverture'] = (df['alignment length'] / df[["query length","subject length"]].max(axis=1)) *100
ndf = df[(df['% identity'] > 70) & (df['evalue'] < 0.01) & (df['% couverture'] >  60)]

df.to_csv('C:/Users/Administrator/Desktop/AMI2B/Gen_Comp/best_hits/best_hits_avec_conditions.csv', sep='\t')


"""
% identity : 70%
evalue : 0,01
condition sur couverture : on divise longueur alignement par longueru max entre querye et subject : 60%


"""


