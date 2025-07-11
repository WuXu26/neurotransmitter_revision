#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import numpy as np
import matplotlib
from matplotlib.collections import LineCollection
from matplotlib import pyplot
from numpy.random import rand
from pylab import figure

from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA

from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import matplotlib.patches as mpatches




class MDS:

    def __init__(self,**kwargs):
        self.setting = kwargs["setting"]
        self.folder = kwargs["outFolder"]
        self.details = kwargs["samplesFile"]
        self.has_columns =  kwargs["has_columns"]

    def assign_colors(self, groups):
        """Randomly assigns colors to the unique groups in the data."""
        colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS).keys()
        color_dict = dict(zip(groups,np.random.choice(colors,len(groups),replace = False)))
        return color_dict

    def get_single_scatter_proxy(self, color, marker):
        """Prepares a proxy plot for data labelling in 3D MDS plots."""
        scatter_proxy = matplotlib.lines.Line2D([0],[0], linestyle="none", c=color, marker = marker)
        return scatter_proxy

    def modify_file(self, file):
        """Replace Feng's output contents of semicolon to comma."""
        text = open(file, "r")
        text = ''.join([i for i in text]) \
                .replace(";", ",")
        f = open(file, "w")
        f.writelines(text)
        f.close()
        df = pd.read_csv(file, index_col=0, header = None)
        groups_protein_dict = dict(zip(self.details['protein'].values, self.details['sampleClass'].values))
        df.index = [groups_protein_dict[i] for i in df.index]
        df.columns = df.index
        return df
        
    def get_mds(self,file, type):
        """Generates MDS in 2D, 3D with and without data labels."""
        if self.has_columns:
            sim_mat = pd.read_csv(file,header = 0,index_col=0)
        else:
            sim_mat = self.modify_file(os.path.join(self.folder,"{}".format(file)))
        #print(self.details)
        samples_file ={}
        for name in self.details['sampleClass'].values:
            samples_file[name.split('-')[0]] = name.split('-')[1]
        
        groups = [str(i.split('-')[0]) for i in sim_mat.columns]
        #Why am i using this file to calculate f? The text is calculated using sim_mat so you that in f, why using sample_details?.
        #groups = [str(i.split('-')[0]) for i in self.details['sampleClass'].values]
        #print('groups: ', groups)
        unique_groups = list(set(groups))
        #print('unique_groups: ', unique_groups)
        color_dict = self.assign_colors(unique_groups)
        #print(color_dict)
        f = [color_dict[str(i)] for i in groups]

        print('--------------------------Start 2-D MDS({})-----------------------------------------'.format(type))
        #Two-dimensional scaling
        mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=123,
                           dissimilarity="precomputed", n_jobs=1)
        pos = mds.fit(sim_mat).embedding_
        fig, ax = plt.subplots()
        ax.scatter(pos[:, 0], pos[:, 1], color=f, s=100, lw=0, label='MDS')
        recs = []
        for i in unique_groups:
            recs.append(mpatches.Rectangle((0,0),1,1,fc=color_dict[i]))
        plt.legend(recs,unique_groups,loc='best', title = 'Protein Groups')
        plt.title('2-D MDS({})'.format(type))
        ax.get_figure().savefig('{}/mds_{}_{}_{}.png'.format(self.folder,'2D', type,self.setting))
        plt.show()
        print('--------------------------End 2-D MDS({})-----------------------------------------'.format(type))


        #print('--------------------------Start 2-D MDS with labels({})-----------------------------------------'.format(type))
        #Two-dimensional scaling
        #mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=123,
        #                   dissimilarity="precomputed", n_jobs=1)
        #pos = mds.fit(sim_mat).embedding_
        #fig, ax = plt.subplots()
        #ax.scatter(pos[:, 0], pos[:, 1], color=f, s=100, lw=0, label='MDS')
        #print(f)
        #for i, txt in enumerate(sim_mat.columns):
        #    print(pos[:, 0][i], pos[:, 1][i], f[i], txt)
        #    ax.annotate(txt, (pos[:, 0][i], pos[:, 1][i]))
            #ax.annotate(txt.split('-')[1], (pos[:, 0][i], pos[:, 1][i]))
        #ax.legend()
        #plt.title('2-D MDS with labels({})'.format(type))
        #ax.get_figure().savefig('{}/mds_{}_{}_{}.png'.format(self.folder,'2D_labels', type,self.setting))
        #plt.show()
        #print('--------------------------End 2-D MDS({})-----------------------------------------'.format(type))


        print('--------------------------Start 2-D MDS({}) without legend-----------------------------------------'.format(type))
        #Two-dimensional scaling
        mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=123,
                           dissimilarity="precomputed", n_jobs=1)
        pos = mds.fit(sim_mat).embedding_
        fig, ax = plt.subplots()
        ax.scatter(pos[:, 0], pos[:, 1], color=f, s=100, lw=0, label='MDS')
        recs = []
        for i in unique_groups:
            recs.append(mpatches.Rectangle((0,0),1,1,fc=color_dict[i]))
        #plt.legend(recs,unique_groups,loc='best', title = 'Protein Groups')
        ax.get_figure().savefig('{}/mds_{}_{}_{}.png'.format(self.folder,'2D_no_legend', type,self.setting))
        #plt.show()
        print('--------------------------End 2-D MDS({}) without legend-----------------------------------------'.format(type))



        print('--------------------------Start 3-D MDS({})-----------------------------------------'.format(type))
        #Three-dimensional scaling
        mds = manifold.MDS(n_components=3, max_iter=3000, eps=1e-9, random_state=123,
                           dissimilarity="precomputed", n_jobs=1)
        pos = mds.fit(sim_mat).embedding_
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(pos[:, 0], pos[:, 1], pos[:,2], c=f, marker='o', s = 60)
        scatter_proxy = [self.get_single_scatter_proxy(color_dict[group], 'o') for group in unique_groups]
        ax.legend(scatter_proxy, \
            unique_groups, \
            scatterpoints = 1, \
            loc = 'lower left', \
            fontsize = 'medium', \
            title = 'Protein Groups', \
            ncol = 4
            )
        plt.title('3-D MDS({})'.format(type))
        ax.get_figure().savefig('{}/mds_{}_{}_{}.png'.format(self.folder,'3D', type,self.setting))
        plt.show()

        print('--------------------------End 3-D MDS({})-----------------------------------------'.format(type))



        print('--------------------------Start 3-D MDS with labels({})-----------------------------------------'.format(type))
        #Three-dimensional scaling with data labels
        m=pos
        fig = figure()
        ax = Axes3D(fig)
        for i in range(len(m)): #plot each point + it's index as text above
            ax.scatter(m[i,0],m[i,1],m[i,2],c=f[i], s = 60) 
            ax.text(m[i,0],m[i,1],m[i,2],  '%s' % (str(sim_mat.columns[i].split('-')[1])), size=14, zorder=25, \
             color='k', va = 'bottom') 
        scatter_proxy = [self.get_single_scatter_proxy(color_dict[group], 'o') for group in unique_groups]        
        ax.legend(scatter_proxy, \
            unique_groups, \
            scatterpoints = 1, \
            loc = 'lower left', \
            fontsize = 'medium', \
            title = 'Protein Groups', \
            ncol = 4
            )
        plt.title('3-D MDS with labels({})'.format(type))
        ax.get_figure().savefig('{}/mds_{}_{}_{}.png'.format(self.folder,'3D_labels', type,self.setting))
        pyplot.show()
        print('--------------------------End 3-D MDS with labels-----------------------------------------'.format(type))


	print('--------------------------Start 3-D MDS({}) without legend-----------------------------------------'.format(type))
        #Three-dimensional scaling
        mds = manifold.MDS(n_components=3, max_iter=3000, eps=1e-9, random_state=123,
                           dissimilarity="precomputed", n_jobs=1)
        pos = mds.fit(sim_mat).embedding_
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(pos[:, 0], pos[:, 1], pos[:,2], c=f, marker='o', s = 60)
        scatter_proxy = [self.get_single_scatter_proxy(color_dict[group], 'o') for group in unique_groups]
        ax.get_figure().savefig('{}/mds_{}_{}_{}.png'.format(self.folder,'3D_no_legend', type,self.setting))
        #plt.show()
        print('--------------------------End 3-D MDS({}) without legend-----------------------------------------'.format(type))


    def visualize_mds(self):
        """Function that is called from outside main.
        Generates MDS in 2D, 3D with and without data labels."""
        # self.get_mds(os.path.join(self.folder, "normal_jaccard_similarity{}.csv".format(self.setting)), 'normal')
        # self.get_mds(os.path.join(self.folder, "generalised_jaccard_similarity{}.csv".format(self.setting)), 'generalised')
        # self.get_mds(os.path.join(self.folder, "wu_jaccard_similarity{}.csv".format(self.setting)), 'wu')
        # self.get_mds(os.path.join(self.folder, "sarika_jaccard1_similarity{}.csv".format(self.setting)), 'sarika')

        #self.get_mds(os.path.join(self.folder, "normal.csv"), 'normal')
        self.get_mds(os.path.join(self.folder, "generalised.csv"), 'generalised')
        #self.get_mds(os.path.join(self.folder, "wu.csv"), 'wu')
        #self.get_mds(os.path.join(self.folder, "sarika.csv"), 'sarika')

        #self.get_mds(os.path.join(self.folder, "cosine.csv"), 'cosine')
