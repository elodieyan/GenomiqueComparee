#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 10:19:52 2020

@author: Chanco
"""

###  Imports

import pandas as pd
import pickle 
import progressbar
from io import StringIO
import sys
import os
from pipeline_functions import *

with open('/Users/Chanco/Desktop/AgroParisTech/3A/AMI2B/Genomique_comparee/'\
          + 'full_pipeline/GenomiqueComparee/list_genomes.pickle','rb') as genomes_file:
    loader = pickle.Unpickler(genomes_file)
    list_genomes = loader.load()

### parameters

analyzed_genomes = list_genomes
Id_perc = 70
ev_treshold = 0.01
cover_perc = 60
small_com = ''

blast_outputs_path = '/Users/Chanco/Desktop/AgroParisTech/3A/AMI2B/Genomique_comparee/full_pipeline/blast_outputs/'
workpath = '/Users/Chanco/Desktop/AgroParisTech/3A/AMI2B/Genomique_comparee/full_pipeline/Analyses/' # chemin où on veut stocker les resultats et fichiers intermediaires

analyse_name = workpath \
            + 'nb_genomes_' + str(len(analyzed_genomes)) \
            + '_id%_' + str(Id_perc)\
            + '_ev_' + str(ev_treshold) \
            + '_cover%_' + str(cover_perc) \
            + small_com + '/'

os.makedirs(analyse_name)

### Main program


# Identify Best Hits

bh_files_path = analyse_name + 'Best_hits/'
os.makedirs(bh_files_path)
for genome1 in analyzed_genomes:
    for genome2 in analyzed_genomes:
        if genome1 == genome2:
            continue
        else:
            file_path = blast_outputs_path + genome1 + '-vs-' + genome2 + '.txt'
            df_best_hit = Best_hit(file_path,Id_perc,ev_treshold,cover_perc)
            new_file_path = bh_files_path + genome1 + '-vs-' + genome2 + '_Best_hit.txt'
            df_best_hit.to_csv(new_file_path, sep='\t', index = False)

# Find BBH
bbh_files_path = analyse_name +'BBH_files/'
os.makedirs(bbh_files_path)
analyzed_genomes2 = analyzed_genomes.copy()



for genome1 in analyzed_genomes:
    for genome2 in analyzed_genomes2:
        if genome1 == genome2:
            continue
        else:
            bbh_file_name = bbh_files_path  + genome1 + '-vs-' + genome2 + '_BBH.csv'
            df_bbh = BBH(bh_files_path,genome1,genome2)
            df_bbh.to_csv(bbh_file_name, sep='\t', index = False)
           
    analyzed_genomes2.remove(genome1)

# Cliques

dico,cliques_nb, genes_min = clique(bbh_files_path)
print(cliques_nb)

# Export results

with open(analyse_name + 'dico.pickle', 'wb') as handle:
    pickle.dump(dico, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
#with open(analyse_name + 'dico_clique.pickle', 'wb') as handle:
#    pickle.dump(dico_clique, handle, protocol=pickle.HIGHEST_PROTOCOL)

summary = 'smaller genome \t nb cliques' + '\n' + str(genes_min) + '\t' +str(cliques_nb)
with open(analyse_name + 'summary.txt', 'w') as file:
    file.write(summary)
                
