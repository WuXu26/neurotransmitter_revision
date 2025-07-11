import csv
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage, cut_tree
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
import seaborn as sns
import os

f=open('generalised.csv','r')
df=pd.read_csv('generalised.csv')
similarity_Matrix=np.empty(shape=(len(df)+1,len(df)+1))
pdbFileindex=[]
lines=f.readlines()
i=0
for line in lines:
    similarity_data=line.split(',')
    pdbFileindex.append(similarity_data[0])
    for j in range(1,len(similarity_data)):
        similarity_Matrix[i][j-1] = similarity_data[j]
        print(j)
        print(i,j)

    i+=1
print(similarity_Matrix)
noOfPdbFiles=len(pdbFileindex)
print(noOfPdbFiles,len(df)+1)
print(pdbFileindex)


#Importing the similarity_Matrix to a csv file
#df_2=pd.DataFrame(similarity_Matrix,columns=pdbFileindex,index=pdbFileindex)
#df_2.to_csv("similarityMatrix_CDRL3.csv",index=True)


similarity_Matrix=abs(similarity_Matrix)
df_3=pd.DataFrame(similarity_Matrix,index=pdbFileindex,columns=pdbFileindex)
print(df_3)
dists = squareform(df_3)
"""
#Separate dataframe to make specific file 'KeyFreq_theta29_dist17.csv'
first_column = list(df.iloc[:, 0])
pdbFileindex_2=[]
for i in range(len(df)):
    pdbFileindex_2.append(pdbFileindex[i]+";"+str(first_column[i]))
df.iloc[:, 0]=pdbFileindex_2
df.to_csv('KeyFreq_theta29_dist17.csv',header=False,index=False)
"""

# Plotting dendrogam (uses average distance)
linkage_matrix = linkage(dists, "average")
dendrogram(linkage_matrix, labels=pdbFileindex)
# plt.xlabel("PDB File(Protein+Drug )", fontsize=15)

plt.ylabel("Distance", fontsize=15)
plt.title("Dendogram", fontsize=15)
plt.savefig('Dendrogram_usr.png')
plt.show()

# Plotting clustermap after keygen code
sns.clustermap(df_3,row_linkage=linkage_matrix, col_linkage=linkage_matrix, annot=False)
plt.savefig('Heatmap_usr.png')
plt.show()

n_clusters = range(1, len(pdbFileindex))
clusters=cut_tree(linkage_matrix,n_clusters )
clusters = np.insert(clusters,clusters.shape[1],range(clusters.shape[0]), axis=1)
clusters=clusters.T

Clustered_Data_File=open("Clustered_data_every_step_usr",'w')
Clustered_Data_File.writelines('Number_of_Cluster----Clustered_PDB_ID\n')

amionAcidEachGroup=open('cluster_usr_Group','w')
amionAcidEachGroup.writelines('Groups\tNo_aa\n')
for row in clusters[::-1]:
    groups={}
    for i,g in enumerate(row):
        if g not in groups:
            groups[g]=set([])
        groups[g].add(pdbFileindex[i])
    clustered_data=list(groups.values())
    Clustered_Data_File.writelines('{}----{}\n'.format(len(clustered_data),clustered_data))
    if len(clustered_data)==2:#group number is here
        groupIndex = 1
        for i in clustered_data:
            amionAcidEachGroup.writelines('Group_{}\n'.format(groupIndex))
            for j in i:
                if j in pdbFileindex:
                    amionAcidEachGroup.writelines('{}\n'.format(j))
                else:
                    print('{} not found'.format(j))


            groupIndex += 1

Clustered_Data_File.close()
amionAcidEachGroup.close()

