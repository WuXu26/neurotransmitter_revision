import glob, os, csv, ntpath,socket,argparse, time, re
import pandas as pd, numpy as np
from collections import Counter
from joblib import Parallel, delayed, cpu_count
from os.path import expanduser
import itertools
import matplotlib as mpl
mpl.use('Agg')
from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
import itertools

__author__ = "Venkata Sarika Kondra"
__version__ = "1.0.1"
__maintainer__ = "Venkata Sarika Kondra"
__email__ = "c00219805@louisiana.edu"

parser = argparse.ArgumentParser(description='Search Keys.')
parser.add_argument('--sample_name', '-sample', metavar='sample_name', \
	default='t1', \
	help='Name of the sample on which this script should be run.')
parser.add_argument('--path', '-path', metavar='path', \
	default=os.path.join(expanduser('~'),'Research', 'Protien_Database', \
		'extracted_new_samples', 'testing'), \
	help='Directory of input sample and other files.')
parser.add_argument('--aacds', '-aacds', metavar='aacds', \
	default='gly,thr,lys',\
	help='Amino acids that are to be searched.')
parser.add_argument('--setting', '-setting', metavar='setting', \
	default='theta29_dist35', \
	help='Setting with theta and distance bin values.')

def get_keys_percents(files, keys, required_percents, group_name):
	#Read keys files only
	print( 'Key Percent Calculation started.')
	common_keys = {}
	common_keys_percents = {}
	start = time.time()
	for file in files:
		#print( ntpath.basename(file)[:4])
		if keys:
			for line in open(file, 'r'):
				keyNo = line.split('\t')[0]
				if keyNo in keys:
					if keyNo in common_keys:
						common_keys[keyNo] += 1
					else:
						common_keys[keyNo] = 1
		else:
			for line in open(file, 'r'):
				keyNo = line.split('\t')[0]
				if keyNo in common_keys:
					common_keys[keyNo] += 1
				else:
					common_keys[keyNo] = 1
	df = pd.DataFrame(common_keys.items(), columns = ['key', '#files_present'])
	df['key_percent'] = (100*df['#files_present']/len(files)).astype(float)	
	for percent in required_percents:
		if percent == 0:
			df_filtered = df
		else:	
			df_filtered = df[df['key_percent'] >= percent]
		common_keys_percents[percent] = df_filtered['key'].values
		#print df_filtered.head(5)
		df_filtered.to_csv(os.path.join(args.path, args.sample_name, \
				args.setting, 'keys_percents_{}_{}.csv'.format(group_name, percent)))
	print( 'Time taken for Key Percent Calculation: ', (time.time() - start)/60)
	return common_keys_percents

def common_keys(files):
	print( 'Common Keys Calculation started.')
	common_keys = []
	start = time.time()
	for file in files:
		print( ntpath.basename(file)[:4])
		keys = []
		for line in open(file, 'r'):
			keys.append(line.split('\t')[0])
		if common_keys:
			common_keys = list(set(common_keys) & set(keys))
		else:
			common_keys = list(set(keys))
		
	print( 'Time taken for Common Keys Calculation: ', (time.time() - start)/60)
	return common_keys

def get_certain_percent_keys_groups():
	samples_file = pd.read_csv(os.path.join(args.path, \
		args.sample_name, 'sample_details.csv'))
	groups = list(set(samples_file['group'].values))

	files = glob.glob(os.path.join(args.path, args.sample_name, args.setting, '*.keys*'))
	samples_file['filename'] = os.path.join(args.path, args.sample_name, args.setting) + '//' +samples_file['protein'].astype(str) + '.keys_'+args.setting
	summary = []
	groups_cKeys = {}
	all_keys_all_groups = {}
	t_all = {}
	print groups
	for group in groups:
		print group
		group_files = samples_file[samples_file['group'] == group]['filename'].values
		group_threshold = list(set(samples_file['percent_threshold'].values))[0]
		keys_percents = get_keys_percents(group_files, None, [group_threshold,0], group)
		groups_cKeys[group] = list(keys_percents[group_threshold])
		all_keys_all_groups[group] = list(keys_percents[0])
		
	for k,v in all_keys_all_groups.items():
		print k, '-', len(v)
	print '------------------------------------------------------------------------------------------'	
	
	df = get_venn_segments(all_keys_all_groups)
	df.to_csv(os.path.join(args.path, args.sample_name, \
					args.setting, 'venn_segment_counts.csv'))
	
	if len(groups) == 3:
		get_venn_3_groups(all_keys_all_groups, \
			'Key distribution among various protease groups\n')


def calculate_intersections(subsets_dict, set_counts, groups):
	only_intersections = {}
	for L in range(len(groups)+1, 0, -1):
	    for subset in itertools.combinations(groups, L):
	    	if subset:
				filtered_set_counts = {k: v for k, v in set_counts.iteritems() if v > len(subset) and set(list(subset)).issubset(k.split('-'))}
				union_list = []
				for grp in filtered_set_counts.keys():
					union_list.extend(subsets_dict[grp])
				only_intersections["-".join(subset)] = len(list(set(subsets_dict["-".join(subset)]) - set(union_list)))
				#only_intersections["-".join(subset)] = list(set(subsets_dict["-".join(subset)]) - set(union_list))
	print only_intersections
	return pd.DataFrame(only_intersections.items(), columns = ['segments', 'counts'])
		
def get_venn_segments(all_keys_all_groups):
	all_combination_intersects_len = {}
	all_combination_intersects = {}
	for L in range(0, len(all_keys_all_groups)+1):
	    for subset in itertools.combinations(all_keys_all_groups, L):
	        common_keys = []
	        for grp in subset:
	        	if len(common_keys) == 0:
	        		common_keys = list(set(all_keys_all_groups[grp]))
	        	else:
	        		common_keys = list(set(common_keys) & set(all_keys_all_groups[grp]))
	        all_combination_intersects["-".join(subset)] = common_keys
	        all_combination_intersects_len["-".join(subset)] = len(subset)
	df = calculate_intersections(all_combination_intersects, all_combination_intersects_len, all_keys_all_groups)
	return df

def get_venn_3_groups(data_dict, title):
	venn3([set(data_dict[data_dict.keys()[0]]), \
		set(data_dict[data_dict.keys()[1]]), \
		set(data_dict[data_dict.keys()[2]])], \
		set_labels = (data_dict.keys()[0], data_dict.keys()[1], data_dict.keys()[2]))
	plt.title(title)
	plt.show()	

if __name__ == '__main__':
	
	start = time.time()
	args = parser.parse_args()
	get_certain_percent_keys_groups()
	#sets = {'A':range(10,20), 'B':range(5,19), 'C':range(15,25), 'D':range(18,30)}
	#get_venn_segments(sets)
	
	
	
