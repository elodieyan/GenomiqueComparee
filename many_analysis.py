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
import random as rd
import sys
import os
import shutil
from pipeline_functions import *



### unvariable parameters

with open('/Users/Chanco/Desktop/AgroParisTech/3A/AMI2B/Genomique_comparee/'\
          + 'full_pipeline/GenomiqueComparee/list_genomes.pickle','rb') as genomes_file:
    loader = pickle.Unpickler(genomes_file)
    list_genomes = loader.load()
blast_outputs_path = '/Users/Chanco/Desktop/AgroParisTech/3A/AMI2B/Genomique_comparee/full_pipeline/blast_outputs/'
workpath = '/Users/Chanco/Desktop/AgroParisTech/3A/AMI2B/Genomique_comparee/full_pipeline/Analyses/' # chemin où on veut stocker les resultats et fichiers intermediaires


def pipeline(genomes_selected,id_perc, ev_tresh, cov_perc, com):
    
    ### variable_parameters
    
    analyzed_genomes = genomes_selected
    Id_perc = id_perc
    ev_treshold = ev_tresh
    cover_perc = cov_perc
    small_com = com
    
    analyse_name = workpath \
            + 'nb_genomes_' + str(len(analyzed_genomes)) \
            + '_id%_' + str(Id_perc)\
            + '_ev_' + str(ev_treshold) \
            + '_cover%_' + str(cover_perc) \
            + small_com + '/'
    
    if os.path.exists(analyse_name):
        shutil.rmtree(analyse_name)
        os.makedirs(analyse_name)
    else:
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
    
    # Export results
    
    with open(analyse_name + 'dico.pickle', 'wb') as handle:
        pickle.dump(dico, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
    #with open(analyse_name + 'dico_clique.pickle', 'wb') as handle:
    #    pickle.dump(dico_clique, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    summary = 'smaller genome \t nb cliques' + '\n' + str(genes_min) + '\t' +str(cliques_nb)
    with open(analyse_name + 'summary.txt', 'w') as file:
        file.write(summary)
                
    
    return (analyse_name, len(genomes_selected),genomes_selected,id_perc, ev_tresh, cov_perc,cliques_nb ,com)


bar = progressbar.ProgressBar(max_value = 6)
i=0
bar.update(i)
for ide in [60,65,70,75,80,85,90]:
    for ev in [0.1,0.05,0.01,0.005,0.001]:
        for cov in [50,60,70,80,90]:
            analyzed_genomes = list_genomes.copy()
            rd.shuffle(analyzed_genomes)
            genomes_selection = analyzed_genomes
            com = ''
            tupl = pipeline(genomes_selection,ide, ev, cov,com)
            df = pd.DataFrame([tupl], columns = ['Analysis_name','nb_genomes','genomes','id %','evalue','cover %','cliques_number','commentary'])
            analyse_table = pd.read_csv('/Users/Chanco/Desktop/AgroParisTech/3A/AMI2B/Genomique_comparee/full_pipeline/Analyses/analyse_table.csv',\
                                            sep = '\t',header=0)
            analyse_table.append(df)
    i += 1
    bar.update(i)
bar.finish()
        
analyse_table.to_csv('/Users/Chanco/Desktop/AgroParisTech/3A/AMI2B/Genomique_comparee/full_pipeline/Analyses/analyse_table.csv',sep='\t',index=False)
    
    





