# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 18:15:55 2020

@author: Administrator
"""



#### imports

import sys
import pandas as pd
from io import StringIO

####

file = open(str(sys.argv[1]),'r')
#bests_hit_file = open(str(sys.argv[2]),'w')

txt = file.read()
querys = txt.split('# BLASTP 2.2.31+')[1:-1]
database_info = querys[0].split('\n')[1:][1]
string =''
string += database_info + '\n'
string += 'query id\tsubject id\t% identity\talignment length\tmismatches\tgap opens\tgaps\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\tquery length\tsubject length' + '\n'
new_name = str(sys.argv[1]).split('/')[-1][:-4] + '_compare.csv'

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

df['evalue'] = df['evalue'].astype(float)
df['% identity'] = df['% identity'].astype(float)
df['% couverture'] = df['% couverture'].astype(float)
ndf = df[['query id','subject id','% identity','evalue','% couverture']]
ndf = ndf[(ndf['% identity'] > 70) & (df['evalue'] < 0.01) & (df['% couverture'] >  60)]

ndf.to_csv(new_name, index = False)
