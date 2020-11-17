# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 10:49:48 2020

@author: Administrator
"""

import pandas as pd

import os

# on indique l'endroit où se trouvent nos fichiers BBH
directory = r'C:\Users\Administrator\Desktop\AMI2B\Gen_Comp\test'

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
                
#dico.to_csv('dico.csv', index = False)   
                
dico_copie = dico
                
# CONCERNANT LES CLIQUES  

# on crée un dictionnaire qui va stocker les cliques
dico_clique = {}

# on crée un dictionnaire pour stocker les génomes déjà parcourus 

genomes_vus = {}
clique_num = 1

# on parcourt les génomes
for genome1 in dico_copie :
    
    for gene1 in dico_copie[genome1] :
        
        if len(dico_copie[genome1][gene1]) == 4 :
            
            for gene2 in dico_copie[genome1][gene1] :
                t = 0
                
                for genome2 in dico_copie :
                    if gene2 in dico_copie[genome2].keys() and dico_copie[genome2][gene2] == dico_copie[genome1][gene1]:
                        t +=1
                    else : 
                        break 
                    i = genome2
                
                #del dico_copie.pop[i][gene2]
                
                
            if t == 4 :
                dico_clique['Clique_'+str(clique_num)] = dico[genome1][gene1]
                clique_num+=1
                        

                            
    
    
#    # pour chaque génome vu, on l'ajoute dans le dictionnaire
#    genomes_vus[genome1]=1
#    
#    # ce qui permet de ne parcourir qu'une fois chaque génome face à un autre
#    # et donc d'éviter les doublons
#    for genome2 in dico.keys() - genomes_vus.keys() :
#        
#        # on parcourt les gènes du génome 1
#        for gene1 in dico[genome1] :
#            if len(dico[genome1][gene1]) == 4 :
#            
#                #et ceux du génome 2
#                for gene2 in dico[genome2] :
#                    
#                    # si les gènes sont les mêmes, on ajoute dans le dico_clique
#                    if dico[genome2][gene2] == dico[genome1][gene1] :
#                        dico_clique[gene2] = dico[genome2][gene2]


#and dico[genome2][gene2] not in dico_clique.values()                   
                    
#import pickle
#                    
#with open('dico.pickle', 'wb') as handle:
#    pickle.dump(dico, handle, protocol=pickle.HIGHEST_PROTOCOL)
#
#with open('dico_clique.pickle', 'wb') as handle:
#    pickle.dump(dico_clique, handle, protocol=pickle.HIGHEST_PROTOCOL)
#    
#                    
#
##file_to_read = open('C:\Users\Administrator\Desktop\AMI2B\Gen_Comp\dico.pickle', 'wb')
#
##dico = pickle.load(file_to_read)
#
#with open(r'C:\Users\Administrator\Desktop\AMI2B\Gen_Comp\dico_clique.pickle','rb') as file:
#    a = pickle.Unpickler(file)
#    dico_clique = a.load()

#print(loaded_dictionary)                    
                    