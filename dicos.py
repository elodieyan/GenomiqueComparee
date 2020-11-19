# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 18:15:07 2020

@author: Elodie
"""
import pandas as pd

import os

# CONSTRUCTION DU DICTIONNAIRE 

# on indique l'endroit où se trouvent nos fichiers BBH
directory = r'C:\Users\Elodie\Documents\Ecole\APT-3A\Master\GenomiqueComparee\Fichiers_bbh'

# on crée le dictionnaire qui contiendra tous les best-its pour chaque gène de chaque génome
dico = {}  

# on parcourt tous les fichiers
for file in os.scandir(directory):
    
    if file.path.endswith(".csv") :
        bbh_file = pd.read_csv(file)
        
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
    
    
        
import pickle
                    
with open('dico.pickle', 'wb') as handle:
    pickle.dump(dico, handle, protocol=pickle.HIGHEST_PROTOCOL)

with open('dico_clique.pickle', 'wb') as handle:
    pickle.dump(dico_clique, handle, protocol=pickle.HIGHEST_PROTOCOL)
