#program to calulate similarity between molecules


import csv
import math
import Bio.PDB
from Bio.PDB import PDBParser
import pandas as pd
import os
from scipy.stats import skew
import numpy as np
from itertools import combinations

df=pd.read_csv('sample_details_remove_incomplete_asp_glu.csv')
PDB_list = df['protein'].to_list()
Chain = df['chain'].to_list()
Drug_name = df['drug_name'].to_list()
Drug_id = df['drug_id'].to_list()
Group = df['group'].to_list()


def calDist(x1, y1, z1, x2, y2, z2):
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)

def Cal_Centroid(res):
    SumX = 0
    SumY = 0
    SumZ = 0
    counter = 0
    for atom1 in residue:
        atomCoord = atom1.get_vector()
        SumX += atomCoord[0]
        SumY += atomCoord[1]
        SumZ += atomCoord[2]
        counter += 1
    Centroid_X_ = SumX / counter
    Centroid_Y_ = SumY / counter
    Centroid_Z_ = SumZ / counter
    return Centroid_X_,Centroid_Y_,Centroid_Z_
def Cal_NA_Centroid(*args):
    Cen_X=args[0]
    Cen_Y=args[1]
    Cen_Z=args[2]
    dist_dict={}
    for atom1 in residue:
        atomCoord = atom1.get_vector()
        dist=calDist(Cen_X, Cen_Y, Cen_Z, atomCoord[0], atomCoord[1], atomCoord[2])
        dist_dict[dist]=f'{atomCoord[0]}_{atomCoord[1]}_{atomCoord[2]}'
    #print(dist_dict)
    #print(min(dist_dict.values()))
    NA_Coord=dist_dict[min(dist_dict.keys())].split("_")
    DA_Coord = dist_dict[max(dist_dict.keys())].split("_")
    NA_Centroid_X_=NA_Coord[0]
    NA_Centroid_Y_ = NA_Coord[1]
    NA_Centroid_Z_ = NA_Coord[2]

    DA_Centroid_X_=DA_Coord[0]
    DA_Centroid_Y_ = DA_Coord[1]
    DA_Centroid_Z_ = DA_Coord[2]

    CM_=np.mean([*dist_dict.keys()])
    CV_=np.var([*dist_dict.keys()])
    CS_=skew([*dist_dict.keys()], axis=0, bias=True)


    return NA_Centroid_X_,NA_Centroid_Y_,NA_Centroid_Z_,DA_Centroid_X_,DA_Centroid_Y_,DA_Centroid_Z_,CM_,CV_,CS_


    #print(dist_dict[min(dist_dict.values())])
def Cal_DPA_Centroid(*args):
    Cen_X=float(args[0])
    Cen_Y=float(args[1])
    Cen_Z=float(args[2])
    dist_dict={}
    for atom1 in residue:
        atomCoord = atom1.get_vector()
        dist=calDist(Cen_X, Cen_Y, Cen_Z, atomCoord[0], atomCoord[1], atomCoord[2])
        dist_dict[dist]=f'{atomCoord[0]}_{atomCoord[1]}_{atomCoord[2]}'
    #print(dist_dict)
    #print(min(dist_dict.values()))
    DPA_Coord = dist_dict[max(dist_dict.keys())].split("_")

    DPA_Centroid_X_=DPA_Coord[0]
    DPA_Centroid_Y_ = DPA_Coord[1]
    DPA_Centroid_Z_ = DPA_Coord[2]

    DAM_=np.mean([*dist_dict.keys()])
    DAV_=np.var([*dist_dict.keys()])
    DAS_=skew([*dist_dict.keys()], axis=0, bias=True)
    return DPA_Centroid_X_,DPA_Centroid_Y_,DPA_Centroid_Z_,DAM_,DAV_,DAS_

def Cal_MVS(*args):
    Cen_X=int(float(args[0]))
    Cen_Y=int(float(args[1]))
    Cen_Z=int(float(args[2]))
    dist_dict={}
    for atom1 in residue:
        atomCoord = atom1.get_vector()
        dist=calDist(Cen_X, Cen_Y, Cen_Z, atomCoord[0], atomCoord[1], atomCoord[2])
        dist_dict[dist]=f'{atomCoord[0]}_{atomCoord[1]}_{atomCoord[2]}'
    #print(dist_dict)
    #print(min(dist_dict.values()))

    M=np.mean([*dist_dict.keys()])
    V=np.var([*dist_dict.keys()])
    S=skew([*dist_dict.keys()], axis=0, bias=True)
    return M,V,S


Dict_drug={}
Output_Folder_Path="/ddnB/work/wxx6941/TSR/code/code/neurotransmitter_revision/output_usr"
OutputFile=open(f'{Output_Folder_Path}/Similarity_drugs.txt','w')
for i in range(len(PDB_list)):
    PDB_ID=PDB_list[i]
    Chain_Name = Chain[i]
    drug_Name = Drug_name[i]
    drug_Id = Drug_id[i]
    group=Group[i]
    Drug=f'{PDB_ID}_{Chain_Name}_{drug_Name}_{drug_Id}_{group}'
    #print(PDB_ID,Chain_Name,drug_Name,drug_Id)
    PDB_File_Path = "/ddnB/work/wxx6941/TSR/code/code/neurotransmitter_revision/PDB_datadir_HRemoved/{}.pdb".format(PDB_ID)
    p = Bio.PDB.PDBParser()
    Structure = p.get_structure('PrimaryStructureChain', PDB_File_Path)
    if drug_Name == 'CL0':
        #drug_Name='CL0'
        drug_Id=f'0{drug_Id}'
    #print(drug_Id)
    #print(drug_Name)

    model = Structure[0]
    for chain in model:
        if chain.id == Chain_Name:

            for residue in chain:

                drug_Identifier = str(residue)[17:18].strip()
                resName = str(residue)[9:12].strip()  # residue Name

                numeric_filter = filter(str.isdigit, str(residue.id))
                Res_Id = "".join(numeric_filter)  # Residue ID


                if drug_Identifier == 'H' and resName == drug_Name and Res_Id == str(drug_Id):
                    Centroid_X,Centroid_Y,Centroid_Z=Cal_Centroid(residue)
                    NA_Centroid_X,NA_Centroid_Y,NA_Centroid_Z,DA_Centroid_X,DA_Centroid_Y,DA_Centroid_Z,CM,CV,CS=Cal_NA_Centroid(Centroid_X,Centroid_Y,Centroid_Z,residue)

                    DPA_Centroid_X,DPA_Centroid_Y,DPA_Centroid_Z,DAM,DAV,DAS=Cal_DPA_Centroid(DA_Centroid_X,DA_Centroid_Y,DA_Centroid_Z,residue)

                    NAM, NAV, NAS = Cal_MVS(NA_Centroid_X,NA_Centroid_Y,NA_Centroid_Z, residue)
                    DPAM, DPAV, DPAS = Cal_MVS(DPA_Centroid_X,DPA_Centroid_Y,DPA_Centroid_Z, residue)
                    Dict_drug[Drug]=f'{CM}_{CV}_{CS}_{NAM}_{NAV}_{NAS}_{DAM}_{DAV}_{DAS}_{DPAM}_{DPAV}_{DPAS}'

print(Dict_drug)
comb_drug = combinations([*Dict_drug], 2)
for i in comb_drug:
    Drug1_moments=Dict_drug[i[0]].split("_")
    Drug2_moments = Dict_drug[i[1]].split("_")
    print(Drug1_moments)
    print(Drug2_moments)
    sum=0
    for k in range(len(Drug1_moments)):
        diff=abs(float(Drug1_moments[k])-float(Drug2_moments[k]))
        sum+=diff
    similarity=1/(1+(sum/12))
    OutputFile.write(f'{i[0]}\t{i[1]}\t{similarity}')
    OutputFile.write('\n')

OutputFile.close()
print('completed successfully')



