# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 18:15:07 2020

@author: Elodie
"""
import pandas as pd

import os

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
        
    
        