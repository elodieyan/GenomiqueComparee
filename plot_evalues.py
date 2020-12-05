# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 17:16:35 2020

@author: Administrator
"""
import pandas as pd
#import pickle 
#import progressbar
from io import StringIO
#import sys
import os
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt
import numpy as np


repertoire = 'C:/Users/Administrator/Desktop/AMI2B/Gen_Comp/Blast_yeasts'

#on prend tous les fichiers csv dans le dossier
fichiers = [f for f in listdir(repertoire) if isfile(join(repertoire, f)) and f.endswith('.csv')] 

#k=0
#liste = []
#for nom in os.listdir(repertoire) :
#    if nom.endswith('.csv'):
#        ofi=open(repertoire+"/"+nom,'r')
#    #le reste, je copie-colle mais tout ce qu'il faut est au dessus
#
#        df_+'k'=pd.read_csv(ofi) 
#        k +=1
#        liste+=[nom]

#on ne veut pas garder les X vs lui même ni les doublons, pas besoin de récuprocité ici        
bons_fichier_igorf =[] 
bons_fichier_cds =[]         
#on va stocker les fichiers deja lus
deja_vu_cds = []
deja_vu_igorf = []
#cette liste va contenir les listes des evalues pour à chaque fois les 15 fichiers
evalues_cds = []
evalues_igorf = []
#une seule longue liste
evalues_cds_ens = []
evalues_igorf_ens = []

#on parcourt les fichiers du dossier
for k in range (0,len(fichiers)) :
    name_1 = fichiers[k]
    #on enlève l'extension du nom
    name_1 = name_1.split('.csv')[0]
    
    #faire truc cds
    #on regarde si c'est un fichier igorf ou cds
    type = name_1.split('_')[3]
    
    #on connait le type on peut enlever le _cds ou _igorf du nom
    name_1 = name_1.split('_cds')[0]
    name_1 = name_1.split('_igorf')[0]
    #on récupère les noms des deux organismes
    b1 = name_1.split('_vs_')[0]
    b2 = name_1.split('_vs_')[1]
    
    #si c'est X contre lui même on ne fait rien 
    if b1 == b2:
        print()
    elif type == 'cds' :
        #on regarde si on a pas deja fait cette paire
        if [b2,b1] in deja_vu_cds:
            print()
        else :
            #si non on ajoute les noms des organismes dans le fichiers deja_vu correspondant au type 
            #et on stocke le nom de ce fichier 
            deja_vu_cds += [[b1,b2]]
            bons_fichier_cds += [fichiers[k]]
            #on ouvre le fichier et on stocke la liste des evalues
            ofi = open(repertoire+"/"+fichiers[k],'r')
            df = pd.read_csv(ofi)
            evalues_cds +=[[df['evalue']]]
            evalues_cds_ens +=[df['evalue']]
    elif type == 'igorf':
        if [b2,b1] in deja_vu_igorf: 
            print()
        else :
            deja_vu_igorf += [[b1,b2]]
            bons_fichier_igorf += [fichiers[k]]
            ofi = open(repertoire+"/"+fichiers[k],'r')
            df = pd.read_csv(ofi)
            evalues_igorf +=[[df['evalue']]]
            evalues_igorf_ens +=[df['evalue']]
 
plt.hist(evalues_cds[2], density=True)    
plt.show()
    
counts,bin_edges = np.histogram(evalues_cds_ens,20,density=True)
bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
plt.plot(bin_centres, counts) 

#myfig, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2)
plt.subplot(2,1,1)
for k in range (0,len(evalues_cds)):
    counts,bin_edges = np.histogram(evalues_cds[k],200,density=True)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    plt.plot(bin_centres, counts) 
    plt.title('CDS VS CDS')
    plt.xlim(0,4000) 
    #plt.show()
plt.subplot(2,1,2)
for k in range (0,len(evalues_igorf)):
    counts,bin_edges = np.histogram(evalues_igorf[k],200,density=True)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    plt.plot(bin_centres, counts) 
    plt.title('IGORF VS IGORF')
    plt.xlim(0,100)       
    plt.show()  
        