# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 16:46:24 2020

@author: Administrator
"""

### Imports

import pandas as pd
#import pickle 
#import progressbar
from io import StringIO
import sys
import os

### 



path_to_file = str(sys.argv[1])
file = open(path_to_file,'r')
txt = file.read()
querys = txt.split('# BLASTP 2.6.0+')[1:-1]
string =''
string += 'query id\tsubject id\t% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\tquery length\tsubject length\tgaps' + '\n'

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
#je ne sais pas pourquoi pour nature subject moi il me fait une liste contenant la liste que l'on veut 
nature_subject = [[len(i.split('_')) for i in list(ndf['subject id'])]][0]

ndf['nature_query'] = nature_query
ndf['nature_subject'] = nature_subject

ndf_cds_vs_cds = ndf[(ndf['nature_query'] == 3) & (ndf['nature_subject'] == 3)]
ndf_igorf_vs_igorf = ndf[(ndf['nature_query'] == 4) & (ndf['nature_subject'] == 4)]             

new_name_cds = path_to_file.split('.txt')[0] + '_cds.csv'        
new_name_igorf = path_to_file.split('.txt')[0] + '_igorf.csv'          
ndf_cds_vs_cds.to_csv(new_name_cds, index = False)
ndf_igorf_vs_igorf.to_csv(new_name_igorf, index = False)
