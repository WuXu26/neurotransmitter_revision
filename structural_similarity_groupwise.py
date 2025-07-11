import os
import argparse
from os.path import expanduser
import pandas as pd
import numpy as np
import collections

parser = argparse.ArgumentParser(description='Hierarchical Classification.')
parser.add_argument('--path', '-path', metavar='path', \
    default=os.path.join(expanduser('~'),'Research', 'Protien_Database', \
        'extracted_new_samples', 'testing'), \
    help='Directory of input sample and other files.')
parser.add_argument('--sample_name', '-sample', metavar='sample_name', \
    default='t1', help='Name of the sample on which this script should be run.')
parser.add_argument('--setting', '-setting', metavar='setting', \
    default='theta29_dist35', \
    help='Name of the setting used for key generation.')
parser.add_argument('--extension', '-extension', metavar='extension', \
    default='_jaccard_similarity_theta29_dist35_NoFeatureSelection_keyCombine0', \
    help='Name of the extension used.')
parser.add_argument('--is_feng', '-is_feng', \
    action='store_true', default=False, \
    help='Enable this option if similarity output is from Fengs code.')


def get_similarity_details(type, writer):
	if args.is_feng:
		sim_mat = pd.read_csv(os.path.join(outFolder, '{}{}.csv'.format(type, args.extension)),\
			index_col = 0)
	else:
		sim_mat = pd.read_csv(os.path.join(outFolder, '{}{}.csv'.format(type, args.extension)),\
		header = 0, index_col = 0)

	cols = pd.isnull(sim_mat).any(1).nonzero()[0]
	sim_mat = sim_mat.drop(sim_mat.columns[cols],axis=1)
	sim_mat = sim_mat.drop(sim_mat.index[cols],axis=0)

	sim_mat = 1 - sim_mat
	upper_tri_arr = np.triu(sim_mat.values).flatten()
	upper_tri_arr = [round(100*elem, 0) for elem in upper_tri_arr]
	
	arr = sim_mat.values[np.nonzero(sim_mat.values)]
	minval = round(100*np.nanmin(arr), 2)	
	mean = round(100*np.nanmean(arr), 2)
	median = round(100*np.nanmedian(arr), 2)
	arr2 = np.where(arr==1, 0, arr)
	maxval = round(100*np.nanmax(arr2), 2)	
	df_sim_freq = pd.DataFrame(collections.Counter(upper_tri_arr).items(), \
		columns = ['similarity', 'freq'])
	df_sim_freq[~df_sim_freq['similarity'].isin([0,100])].to_excel(writer,sheet_name=type)	
	return((minval, maxval, mean, median, type), collections.Counter(upper_tri_arr))

if __name__ == '__main__':
	args = parser.parse_args()
	outFolder = os.path.join(args.path, args.sample_name, args.setting)
	sim_types = ['generalised', 'normal', 'wu', 'sarika']
	writer = pd.ExcelWriter(os.path.join(outFolder, 'similarity_frequency.xlsx'), \
		engine='xlsxwriter')
	sim_details =[]
	all_frequencies = []
	for sim_type in sim_types:
		res = get_similarity_details(sim_type, writer)
		sim_details.append(res[0])
		res[1]['name'] = sim_type
		all_frequencies.append(res[1])
	#sim_details = [get_similarity_details(sim_type, writer) for sim_type in sim_types]
	print(pd.DataFrame(sim_details, \
		columns = ['minimum', 'maximum', 'mean', 'median', 'type']))
	pd.DataFrame(sim_details, \
		columns = ['minimum', 'maximum', 'mean', 'median', 'type']).to_csv(
		os.path.join(outFolder, 'similarity_details.csv'))
	df_sim_freq = pd.DataFrame(all_frequencies)
	df_sim_freq.drop([0.0,100.0], axis = 1).to_excel(writer,sheet_name='all')
	
