
import pandas as pd
import numpy as np
import csv
Output_Folder_Path="/ddnB/work/wxx6941/TSR/code/code/neurotransmitter_revision/output_usr"
file=open(f'{Output_Folder_Path}/Similarity_drugs.txt','r')
lines=file.readlines()
dict={}

for line in lines:
    PDB=line.strip().split()[0]
    similarity=line.strip().split()[2]
    PDB_=line.strip().split()[1]
    if PDB not in dict:
        dict[PDB]=[]
        dict[PDB].append(similarity)

    else:
        dict[PDB].append(similarity)

similarity_matrix=np.zeros((len(dict)+1,len(dict)+1))
#similarity_matrix=similarity_matrix.astype(int)
for PDB in dict:
    print(PDB)
print(len(similarity_matrix))
i=0
k=1
for PDB in dict:
    j=k
    for data in dict[PDB]:
        print(i,j)
        print(j,i)
        similarity_matrix[i][j]= 1-float(data)
        similarity_matrix[j][i] = 1-float(data)


        j+=1
    i+=1
    k+=1
PDB_list=list(dict.keys())
PDB_list.append(PDB_)

DataFrame_Index=[]
for i in range(len(PDB_list)):
    DataFrame_Index.append(PDB_list[i].upper())


df=pd.DataFrame(similarity_matrix,index=DataFrame_Index)
print(df)
#first_column = list(df.iloc[:, 0])

df.to_csv(f"{Output_Folder_Path}/generalised.csv", header=False, index=DataFrame_Index)

print('completed_successfully')



