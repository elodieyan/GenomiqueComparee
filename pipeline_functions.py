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
    df_BBH = pd.DataFrame(columns = ['query id','subject id','% identity','evalue','% couverture'])
    for i in range (0,len(file1.index)):
        if (file2.loc[file2['subject id'] == file1.loc[file1.index[i],'query id']]).empty == False:
            x = (file2.loc[file2['subject id'] == file1.loc[file1.index[i],'query id']])
            for k in range (0,len(x.index)):
                if x.loc[x.index[k],'query id'] == file1.loc[file1.index[i],'subject id']:
                    df_BBH.loc[len(df_BBH.index)] = file1.loc[i,:].values
    return df_BBH
    
    

def clique(repertory_path):
        # on crée le dictionnaire qui contiendra tous les best-its pour chaque gène de chaque génome
    dico = {}  
    
    # on parcourt tous les fichiers
    for file in os.scandir(repertory_path + 'BBH_files/'):
        
        if file.path.endswith(".csv") :
            bbh_file = pd.read_csv(file,header=0,sep='\t')
            
            # pour chaque fichier, on vérifie si le génome est déjà une du dictionnaire
            # si ce n'est pas le cas, on crée la clé qui sera du nom du génome (on constitue en fait un sous-dictionnaire)
            if bbh_file['query id'][0].split('_')[0] not in dico:
                dico[bbh_file['query id'][0].split('_')[0]] = {}  
            
            #on parcourt le fichier
            for i in range (0,len(bbh_file.index)):
                
                # pour chaque gène, on vérifie s'il est déjà une clé du dictionnaire
                # si c'est le cas, on ajoute le gène best-hit de l'autre génome comme valeur pour cette clé
                if bbh_file['query id'][i] in dico[bbh_file['query id'][0].split('_')[0]]:
                    dico[bbh_file['query id'][0].split('_')[0]][bbh_file['query id'][i]].append(bbh_file['subject id'][i])
                
                # si ce n'est pas le cas, on crée la clé puis on ajoute le gène best-hit de l'autre génome comme valeur pour cette clé
                elif bbh_file['query id'][i] not in dico[bbh_file['query id'][0].split('_')[0]] :
                    dico[bbh_file['query id'][0].split('_')[0]][bbh_file['query id'][i]] = [] 
                    dico[bbh_file['query id'][0].split('_')[0]][bbh_file['query id'][i]].append(bbh_file['subject id'][i])
        # CONCERNANT LES CLIQUES 
    # recherche du genome avec le moins de gènes, 
    # cela servira de référence pour la recherche des cliques
    genes_min = 10000
    for genome in dico :
        if len(dico[genome]) < 10000 :
            genes_min = len(dico[genome])
            genome_min = genome
    print(genome_min, len(dico[genome_min]))                
    
    # on crée un dictionnaire qui va stocker les cliques
    dico_clique = {}
    clique_num = 1
    
    for gene1 in dico[genome_min]:
        if len(dico[genome_min][gene1]) == len(dico) :
            for genome in dico.keys()-genome_min:
                for gene2 in dico[genome]:
                    if len(dico[genome][gene2]) == len(dico) and dico[genome][gene2] == dico[genome_min][gene1] and dico[genome][gene2] not in dico_clique.values():
                       dico_clique['Clique_'+str(clique_num)] = dico[genome][gene2]
                       clique_num+=1 
                    else:
                        continue
        else:
            continue
    
    return dico,dico_clique