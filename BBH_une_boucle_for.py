# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 15:30:11 2020

@author: Administrator
"""

import sys
import pandas as pd
from io import StringIO



#faudra voir comment on séléctionne les deux bons fichiers automatiquement
#name_1 = 'Escherichia_coli_536-vs-Escherichia_coli_55989.csv'
#name_2 = 'Escherichia_coli_55989-vs-Escherichia_coli_536.csv'
#file_1 = pd.read_csv(name_1) #colonne 1 : 536, colonne 2 : 55989
#file_2 = pd.read_csv(name_2) #colonne 1 : 55989, colonne 2 : 536

name_1 = str(sys.argv[1])
#name_1 = 'Escherichia_coli_536-vs-Escherichia_coli_55989_compare.csv'
file_1 = pd.read_csv(name_1)
name_1 = name_1.split('_compare.csv')[0]

b1 = name_1.split('-vs-')[0]
b2 = name_1.split('-vs-')[1]

name_2 = b2 + '-vs-' + b1 +'_compare.csv'
file_2 = pd.read_csv(name_2)

#faudra changer les noms des colonnes en fonction des noms des bactérie, ici 536 puis 55989
#on est obligé de mettre le même nom de colonne que pour les file sinon on peut pas copier coller les lignes
file_BBH = pd.DataFrame(columns = ['query id','subject id','% identity','evalue','% couverture'])

#on veut vérif
for i in range (0,len(file_1.index)):
    #on parcourt les lignes du premier fichier,
    #on va regarder pour chaque ligne quel gène est dans la deuxième colonne 
    #puis on va voir dans le deuxième fichier si pour ce gène la dans la permière colonne 
    #le gène dans la deuxième est celui dans la première du fichier 1
    
    
    #on regarde ou est ce que dans le deuxième fichier le subject est le même que celui pour le query
    if (file_2.loc[file_2['subject id'] == file_1.loc[file_1.index[i],'query id']]).empty == False:
        #x peut contenir plusieurs lignes
        x = (file_2.loc[file_2['subject id'] == file_1.loc[file_1.index[i],'query id']])
        #on parcourt les lignes de x 
        for k in range (0,len(x.index)):
            #on regarde si il y a réciprocité
            if x.loc[x.index[k],'query id'] == file_1.loc[file_1.index[i],'subject id']:
                #si oui on ajoute cette ligne dan s le fichier contenant les bbh 
                file_BBH.loc[len(file_BBH.index)] = file_1.loc[i,:].values
                
new_name = name_1 + '_BBH.csv'                
file_BBH.to_csv(new_name, sep='\t', index = False)

        

    
