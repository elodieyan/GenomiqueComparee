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

####

file = open(str(sys.argv[1]),'r')
bests_hit_file = open(str(sys.argv[2]),'w')
file = open('Blast_output/Escherichia_coli_536-vs-Escherichia_coli_55989.txt','r')
bests_hit_file = open('Escherichia_coli_536-vs-Escherichia_coli_55989_compare.txt','w')
txt = file.read()
querys = txt.split('# BLASTP 2.11.0+')[1:-1]
database_info = querys[0].split('\n')[1:][1]
bests_hit_file.write(database_info + '\n')

for q in querys:
    lines = q.split('\n')[1:]
    if lines[2] == '# 0 hits found':
        #best_hit = lines[0][9:]+ '\t' +'no hits found'
        #querys_with_no_hits_found.append(lines[0][9:])
        continue
    else:
        best_hit = lines[4]
        bests_hit_file.write(best_hit + '\n')
    
### Critères de séléction

"""
% identity : 70%
evalue : 0,01
à définir seuil sur alignment length


"""