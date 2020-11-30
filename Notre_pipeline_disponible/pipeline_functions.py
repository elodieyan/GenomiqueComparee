#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 10:17:23 2020

@author: Chanco
"""

###  Modules

import pandas as pd
import pickle 
import progressbar
from io import StringIO
import sys
import os

### Functions
    
def Best_hit(path_to_file,Id_perc, ev_treshold, cover_perc):
    
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
    df['% couverture'] = (df['alignment length'] / df[["query length","subject length"]].max(axis=1)) *100
    df['evalue'] = df['evalue'].astype(float)
    df['% identity'] = df['% identity'].astype(float)
    df['% couverture'] = df['% couverture'].astype(float)
    ndf = df[['query id','subject id','% identity','evalue','% couverture']]
    ndf = ndf[(ndf['% identity'] > Id_perc) & (df['evalue'] < ev_treshold) & (df['% couverture'] >  cover_perc)]
    
    return ndf
    
    

def BBH(Best_hit_path,genome1,genome2):
    file1 = pd.read_csv(Best_hit_path + genome1 + '-vs-' + genome2 + '_Best_hit.txt', header=0, sep='\t')
    file2 = pd.read_csv(Best_hit_path + genome2 + '-vs-' + genome1 + '_Best_hit.txt', header=0, sep='\t')
    dic1_keys = list(file1['query id'])
    dic1_values = list(file1['subject id'])
    
    dic2_keys = list(file2['query id'])
    dic2_values = list(file2['subject id'])
    
    dic1 = dict(zip(dic1_keys,dic1_values))
    dic2 = dict(zip(dic2_keys,dic2_values))
    querys = []
    subjects = []
    t = 0
    for key in dic1_keys:
        t+=1
        if dic1[key] in dic2_keys and dic2[dic1[key]] == key:
            querys.append(key)
            subjects.append(dic1[key])
        else:
            continue
        
    df_BBH = pd.DataFrame({'query id':querys,'subject id':subjects})
    return df_BBH
        
    

def clique(repertory_path):
    dico = {}

    for file in os.scandir(repertory_path):
        
        filename = file.name
        genome1_name = filename.split('-vs-')[0]
        genome2_name = filename.split('-vs-')[1].split('_BBH.csv')[0]
        bbh_file = pd.read_csv(file, header=0,sep='\t')
        
        for genome in [genome1_name,genome2_name]:
            if genome not in dico:
                dico[genome] = {}
        
        for gene_query,gene_subject in zip(list(bbh_file['query id']),list(bbh_file['query id'])):
            if gene_query not in dico[genome1_name]:
                dico[genome1_name][gene_query] = [gene_subject]
            else:
                dico[genome1_name][gene_query].append(gene_subject)
                
            if gene_subject not in dico[genome2_name]:
                dico[genome2_name][gene_subject] = [gene_query]
            else:
                dico[genome2_name][gene_subject].append(gene_query)
        
    genes_min = 10000
    for genome in dico :
        if len(dico[genome]) < genes_min :
            genes_min = len(dico[genome])
            genome_min = genome   
    
    # on crée un dictionnaire qui va stocker les cliques
    cliques_nb = 0
    dic = dico[genome_min]
    l = list(dic.values())
    for i in l:
        if len(i) == len(dico) - 1:
            cliques_nb += 1
    
    return dico, cliques_nb, genes_min
