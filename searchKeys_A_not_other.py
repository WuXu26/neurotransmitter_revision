#!/usr/bin/env python
#searchKeys.py
"""
	Searches for specific keys or amino acids or positions and retrieves 
	their details from triplets files
"""

import glob, os, csv, ntpath,socket,argparse, time, re
import pandas as pd, numpy as np
from collections import Counter
from joblib import Parallel, delayed, cpu_count
from os.path import expanduser
import itertools
import matplotlib as plt

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
parser.add_argument('--keys', '-keys', metavar='keys', \
	default='6362130',\
	help='Key that is to be searched.')
parser.add_argument('--keys_file', '-keys_file', metavar='keys_file', \
	default=None,\
	help='Path to the file where the keys to be searched are are located')
parser.add_argument('--aacds', '-aacds', metavar='aacds', \
	default='gly,thr,lys',\
	help='Amino acids that are to be searched.')
parser.add_argument('--setting', '-setting', metavar='setting', \
	default='theta29_dist35', \
	help='Setting with theta and distance bin values.')
parser.add_argument('--exclude', '-exclude', metavar='exclude', \
	default=True, \
	help='Exclude certain keys from when extracting triplets details.')
parser.add_argument('--search_mode', '-search_mode', metavar='search_mode', \
	default=0, \
	help='0 if Key search, 2 if amino acid search and 3 if sequence identification.')
parser.add_argument('--identifying_pattern', '-identifying_pattern', \
	metavar='identifying_pattern', \
	default='leu1', \
	help='pllop if ploop, leu1 if lxxll and leu2 if llxxl.')
parser.add_argument('--low_freq', '-low_freq', metavar='low_freq', \
	default=5000, \
	help='Keys with frequencies less than this number are considered.')
parser.add_argument('--dist_less', '-dist_less', metavar='dist_less', \
	default=5000, \
	help='Keys with maxDist less than this number are considered.')
parser.add_argument('--use_common_keys', help='Use Common Keys',
    action='store_true')
parser.add_argument('--no_of_xs', help='Number of Xs in general pattern mode',
    default=3)
parser.add_argument('--gen_pattern_aacds', help='amino acids in general pattern mode',
    default='pro,thr')
parser.add_argument('--has_pattern', help='Check if the sample_details file has pattern',
    action='store_true')
parser.add_argument('--a_percent', '-a_percent', metavar='a_percent', \
	default=90, \
	help='Keys with  document frequencies less than this number are considered.')
parser.add_argument('--b_percent', '-b_percent', metavar='b_percent', \
	default=10, \
	help='Keys with frequencies less than this number are considered.')
parser.add_argument('--group_column_name', '-group_column_name', metavar='group_column_name', \
    default='group', help='Select a column name to make it the Y value in classifier.')
parser.add_argument('--percent_column_name', '-percent_column_name', metavar='percent_column_name', \
    default='percent_threshold', help='column name containing pdb ids. protein or PDB ID.')

def search_by_key(fileName, file, keys, iskeys_file, outFolder):
	print(fileName)
	search_records = []
	if iskeys_file == 1:
		line_count = len(file)
		df = file[file['key'].isin(keys)]

	else:
		line_count = 0
		for line in file:
			line_count += 1
			line_splitted=line.split('\t')
			if line_splitted[0].strip() in keys:
				search_records.append(line_splitted)	
		df = pd.DataFrame(search_records,columns = column_names)
	df.to_csv(os.path.join(outFolder, \
				fileName +'_key_search_common_keys.csv'))		
	return (fileName, line_count, len(set(df['key'].values)), len(df['key'].values))


def search_aacd(file,aacds, req_pos):
	search_records = []	
	pattern_line = []
	for line in file:
		line_aacds = {}
		line_splitted=line.split('\t')
		line_aacd_arr = [line_splitted[1].strip().upper(), \
			line_splitted[3].strip().upper(), line_splitted[5].strip().upper()]	
		if sorted(line_aacd_arr) == sorted(aacds):
			search_records.append(line_splitted)
			if req_pos:
				line_pos_arr = [line_splitted[2].strip().upper(), \
					line_splitted[4].strip().upper(), \
					line_splitted[6].strip().upper()]
				if all(elem in req_pos for elem in line_pos_arr):
					pattern_line.append(line)
			
	return pd.DataFrame(search_records,columns = column_names), pattern_line

def identify_pattern_by_regex(files,column_names, patterns):
	matched = []
	all_pattern = generate_sequences_from_triplets_by_FR(file)
	for pattern in patterns:
		print(re.findall(pattern, all_pattern))
		matched.append(pattern)
	return ntpath.basename(file)[:4], matched

def identify_pattern_by_normal(file, pattern,count_of_x, group):
	print(ntpath.basename(file)[:4])
	matched = []
	# If file size is greater than 5 GB pandas is throwing Memory error
	# But if you always read line-by-line it is very slow. 
	# Hence using two way reading
	if float(os.path.getsize(file)/ 2**30) > 4 :
		print('File Size too big')
		return ntpath.basename(file)[:4], ['File Size too big']
		#all_pattern = generate_sequences_from_triplets_by_FR(file)	
	#all_pattern = generate_sequences_from_triplets_by_DF(file)
	all_pattern = generate_sequences_from_triplets_by_FR(file)	
	print all_pattern
	if all_pattern == '':
		return (ntpath.basename(file)[:4], ['Invalid Literal in file.'])
	i = 0
	while i < len(all_pattern):
		if (all_pattern[i].split('_')[0] == pattern[0]) and ( i + count_of_x + len(pattern) < len(all_pattern)):
			#print('checking this: {}'.format(all_pattern[i:i + count_of_x + len(pattern)]))
			try:
				if group == 'ploop':
					if (all_pattern[i + count_of_x +1].split('_')[0] == pattern[1]) \
						and (all_pattern[i+count_of_x +2].split('_')[0] == pattern[2]) \
						and (all_pattern[i+count_of_x +3].split('_')[0] == pattern[3]):
						print( 'found', all_pattern[i:i + count_of_x + len(pattern)])
						matched.append(all_pattern[i:i + count_of_x + len(pattern)])
				if group == 'leu1':
					if (all_pattern[i + count_of_x +1].split('_')[0] == pattern[1]) \
						and (all_pattern[i+count_of_x +2].split('_')[0] == pattern[2]):
						print( 'found', all_pattern[i:i + count_of_x + len(pattern)])
						matched.append(all_pattern[i:i + count_of_x + len(pattern)])
				if group == 'leu2':
					if (all_pattern[i + 1].split('_')[0] == pattern[1]) \
						and (all_pattern[i+count_of_x +2].split('_')[0] == pattern[2]):
						print( 'found', all_pattern[i:i + count_of_x + len(pattern)])
						matched.append(all_pattern[i:i + count_of_x + len(pattern)])
				if group == 'gen':
					if (all_pattern[i + count_of_x + 1].split('_')[0] == pattern[1]) :
						print( 'found', all_pattern[i:i + count_of_x + len(pattern)])
						matched.append(all_pattern[i:i + count_of_x + len(pattern)])
			except ValueError:
				return (ntpath.basename(file)[:4], ['Invalid Literal in file.'])

		i += 1
	return (ntpath.basename(file)[:4], matched)

def generate_sequences_from_triplets_by_DF(file):
	print(ntpath.basename(file)[:4])
	all_pattern = []
	df = pd.read_table(file, delimiter = '\t', names = column_names)
	a_0 = pd.Series(df.aa0.values,index=df.pos0).to_dict()
	a_1 = pd.Series(df.aa1.values,index=df.pos1).to_dict()
	a_2 = pd.Series(df.aa2.values,index=df.pos2).to_dict()
	a_0.update(a_1) 
	a_0.update(a_2) 
	pos_aa_dict = sorted(a_0.items())
	for key, value in pos_aa_dict:
		all_pattern .append(str(value) + '_' + str(key))
	return all_pattern

def generate_sequences_from_triplets_by_FR(file):
	print( ntpath.basename(file)[:4])
	all_pattern = []
	pos_aa = {}
	for line in open(file,'r'):
		line_splitted=line.split('\t')
		try:
			if line_splitted[2].strip() not in pos_aa.keys():
				pos_aa[int(line_splitted[2].strip())] = line_splitted[1].strip().upper()
			if line_splitted[4].strip() not in pos_aa.keys():
				pos_aa[int(line_splitted[4].strip())] = line_splitted[3].strip().upper()
			if line_splitted[6].strip() not in pos_aa.keys():
				pos_aa[int(line_splitted[6].strip())] = line_splitted[5].strip().upper()
		except ValueError:
				return ''
	for key, value in sorted(pos_aa.items()):
		all_pattern .append(str(value) + '_' + str(key))
	return all_pattern

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

def get_keys_percents(files, keys, required_percents, group_name):
	"""For the specified group calcualte all the common keys,
	keys and their respective percents in the group.
	Parametres:
	files:	List of Keys file paths
	keys:	If any specific list of keys to be verified.
	required_percents: [group_hreshold, zero_threshold to allow all keys]
	group_percent:	 Minimum allowable percent of keys in this group.
	other_groups_percent: Maximun allowable percent of keys in this 
		group when this group is part of oneVs all, in the all groups part.
	group_name:	The name of the current group being analysed.

	Returns:
	common_keys_percents: A dictionary with the percentage of one Vs all groups.
			Specific minimum percent required is specified in sample_details file
			Not group percent is mentioned in script argument for 'b_percent'.
	Saves this group both keys- their percents in seperate files
	"""
	#Read keys files only
	# print( 'Key Percent Calculation started. Group Name: {}, #files in group: {} '.\
	# 	format(group_name, len(files)))
	common_keys = {}
	common_keys_percents = {}
	start = time.time()
	for file in files:
		#print( ntpath.basename(file))
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
	df['key_percent'] = (100*df['#files_present']).div(len(files))#.astype(float)	
	df['total_files'] = len(files)
	#plot_group_analysis(df)

	for percent in required_percents:
		if percent == 0:
			#this is to allow all keys without filtering
			df_filtered = df
		else:	
			df_filtered = df[df['key_percent'] >= percent]
		common_keys_percents[percent] = df_filtered['key'].values
		#print df_filtered.head(5)
		df_filtered.to_csv(os.path.join(args.path, args.sample_name, \
				args.setting, 'keys_percents_{}_{}.csv'.format(group_name, percent)))
		print('File saved at {}'.format(os.path.join(args.path, args.sample_name, \
				args.setting, 'keys_percents_{}_{}.csv'.format(group_name, percent))))
	#print( 'Time taken for Key Percent Calculation: ', (time.time() - start)/60)
	return common_keys_percents

def calculate_low_less15_freq_from_common_keys(outFolder, req_low_freq, req_max_dist):	
	summary = []
	commom_key_freqs = Counter()
	if int(args.search_mode) == 5:
		if not os.path.exists(os.path.join(outFolder,'dist_freq_keys')):
			os.makedirs(os.path.join(outFolder,'dist_freq_keys'))
		writer = pd.ExcelWriter(os.path.join(outFolder, 'dist_freq_keys',\
		'all_common_key_distribution.xlsx'), \
		engine='xlsxwriter')
		writer_low_freq = pd.ExcelWriter(os.path.join(outFolder, 'dist_freq_keys',\
			'only_low_freq{}_key_distribution.xlsx'.format(req_low_freq)), \
			engine='xlsxwriter')
		writer_dist_less15 = pd.ExcelWriter(os.path.join(outFolder, 'dist_freq_keys',\
			'freq_less{}_dist_less{}_triplets.xlsx'.\
			format(str(req_low_freq), str(req_max_dist))), engine='xlsxwriter')
		common_keys_files = glob.glob(os.path.join(outFolder, \
		'*.triplets*'))
	else:
		writer = pd.ExcelWriter(os.path.join(outFolder,\
		'all_common_key_distribution.xlsx'), \
		engine='xlsxwriter')
		writer_low_freq = pd.ExcelWriter(os.path.join(outFolder,\
			'only_low_freq{}_key_distribution.xlsx'.format(req_low_freq)), \
			engine='xlsxwriter')
		writer_dist_less15 = pd.ExcelWriter(os.path.join(outFolder,\
			'freq_less{}_dist_less{}_triplets.xlsx'.\
			format(str(req_low_freq), str(req_max_dist))), engine='xlsxwriter')
		common_keys_files = glob.glob(os.path.join(outFolder, \
			'*_key_search_common_keys.csv'))

	for file in common_keys_files:
		print( ntpath.basename(file)[:4])
		if int(args.search_mode) == 5:
			df = pd.read_csv(file, names = column_names, delimiter = '\t')
		else:
			df = pd.read_csv(file, delimiter = ',')
		
		x = Counter(df['key'].values)
		ckeys_freqs_df = pd.DataFrame(sorted(x.items(), key=lambda pair: pair[1], \
			 reverse=True), columns = ['key', 'freq'])
		ckeys_freqs_df.to_excel(writer,sheet_name=ntpath.basename(file)[:4])
		#Low Frequency Common Keys
		low_freqs = ckeys_freqs_df[ckeys_freqs_df['freq'] <= int(req_low_freq)]
		low_freqs.to_excel(writer_low_freq,sheet_name=ntpath.basename(file)[:4])
		df_low_freqs_triplets = df[ df['key'].isin(set(low_freqs['key'].values))]
		#Distance less than 15
		df_distance_less15 = df_low_freqs_triplets.loc[df_low_freqs_triplets['distance'] <= float(req_max_dist)]
		#print df_distance_less15
		if int(args.search_mode) == 5:
			df_distance_less15.to_csv((os.path.join(outFolder,'dist_freq_keys', '{}_freq_less{}_dist_less{}.csv'.\
				format(ntpath.basename(file)[:4],str(req_low_freq), str(req_max_dist)))))
		else:
			df_distance_less15.to_excel(writer_dist_less15,sheet_name=ntpath.basename(file)[:4])
		summary.append((ntpath.basename(file)[:4], len(set(low_freqs['key'].values)),len(low_freqs), \
			len(set(df_distance_less15['key'].values)), len(df_distance_less15)))
		print (ntpath.basename(file)[:4], len(set(low_freqs['key'].values)),len(low_freqs), \
			len(set(df_distance_less15['key'].values)), len(df_distance_less15['key']))
		print (ntpath.basename(file)[:4], len(set(low_freqs['key'].values)),len(low_freqs), \
			len(set(df_low_freqs_triplets['key'].values)), len(df_low_freqs_triplets['key']))
		#Overall Common Keys
		if commom_key_freqs:
			commom_key_freqs = commom_key_freqs + x
		else:
			commom_key_freqs = x
	if int(args.search_mode) != 5:
		pd.DataFrame(sorted(commom_key_freqs.items(), key=lambda pair: pair[1], reverse=True)). \
			to_excel(writer,sheet_name='Summary')

		df_merge = pd.read_csv(os.path.join(outFolder, \
			'summary_common_keys.csv'), dtype = 'unicode').merge(pd.DataFrame(summary, \
				columns = ['fileName', '# of low freq common keys', \
				'# of low freq common keys(with freq)',\
				'# of low freq common keys with maxdist <= {}'.format(str(req_max_dist)), \
				'# of low freq common keys with maxdist <= {}(with freq)'.format(str(req_max_dist))]),\
			on = 'fileName')
		df_merge.to_csv(os.path.join(outFolder, 'summary_common_keys_merged.csv'))
	writer.close()
	writer_low_freq.close()
	writer_dist_less15.close()

def get_common_keys_groups(group_column_name):
	writer = pd.ExcelWriter(os.path.join(args.path, args.sample_name, \
					args.setting, 'common_keys_group_specific.xlsx'), \
					engine = 'xlsxwriter')
	samples_file = pd.read_csv(os.path.join(args.path, \
		args.sample_name, 'sample_details.csv'))
	groups = list(set(samples_file[group_column_name].values))
	files = glob.glob(os.path.join(args.path, args.sample_name, args.setting, '*.keys*'))
	# all_common = common_keys(files)
	# pd.DataFrame(all_common, columns = ['key'])\
	# 					.to_excel(writer,sheet_name='all_common_keys')
	samples_file['filename'] = os.path.join(args.path, args.sample_name, args.setting) + '//' +samples_file['protein']+ '.keys_' + args.setting
	summary = []
	groups_cKeys = {}
	t_all = []
	for group in groups:
		group_files = samples_file[samples_file[group_column_name] == group]['filename'].values
		groups_cKeys[group] = common_keys(group_files)
	print(len(groups_cKeys))	
	for group in groups:
		only_group_keys = groups_cKeys[group]
		print 'All keys of group ', group, '-',len(only_group_keys)
		t = []
		for group2 in groups:
			if group != group2:
				print group, '-', group2
				only_group_keys = list(set(only_group_keys) - set(groups_cKeys[group2]))
				print(len(only_group_keys))
				t.append(len(only_group_keys))
			else:
				t.append(len(groups_cKeys[group]))
		t_all.append(t)

		pd.DataFrame(list(only_group_keys), \
						columns = ['key'])\
						.to_excel(writer,sheet_name=group)
		summary.append((group, len(groups_cKeys[group]), len(only_group_keys)))
	print t_all
	pd.DataFrame(summary, columns = [group_column_name, 'All keys in group', 'Only group keys without common'])\
						.to_excel(writer,sheet_name='summary')
	writer.close()

def plot_group_analysis(df):
	df_percent_grouping = df[['key_percent']].groupby(['key_percent']).agg(['count'])
	print(df_percent_grouping)

def get_certain_percent_keys_groups(percent_column_name, group_column_name):
	writer = pd.ExcelWriter(os.path.join(args.path, args.sample_name, \
					args.setting, 'common_keys_group_specific.xlsx'), \
					engine = 'xlsxwriter')
	samples_file = pd.read_csv(os.path.join(args.path, \
					args.sample_name, 'sample_details.csv'))
	#Get the groups that are to be studied
	groups = list(set(samples_file[group_column_name].values))
	print('Groups in analysis: {}'.format(groups))
	#Use the keys files for further analysis
	files = glob.glob(os.path.join(args.path, args.sample_name, args.setting, '*.keys*'))
	#Add full path of the keys file in the sample details file object
	samples_file['filename'] = os.path.join(args.path, args.sample_name, args.setting) + '//' +samples_file['protein'] + '.keys_'+args.setting
	summary = []
	inter_summary = []
	groups_cKeys = {}#Group specific keys above/equal to group threshold
	all_keys_all_groups = {}#Store all keys of all groups
	t_all = {}
	for group in groups:
		print ("Current Group in Analysis: ", group)
		#Get all the files belonging to this group
		group_files_df = samples_file[samples_file[group_column_name] == group]
		group_files = group_files_df['filename'].values
		#Read the minimum allowable group threshold from the sample details file
		group_threshold = list(set(group_files_df[percent_column_name].values))[0]
		#Passing 0 to get all keys of the group
		keys_percents = get_keys_percents(group_files, None, [group_threshold,0], group)
		groups_cKeys[group] = list(keys_percents[group_threshold])
		all_keys_all_groups[group] = list(keys_percents[0])
		inter_summary.append((group, group_threshold, 
			len(groups_cKeys[group]), len(all_keys_all_groups[group])))
	print("-"*120)
	print("Keys in all groups according to their group thresholds.")
	print(pd.DataFrame(inter_summary, 
		columns = ['Group', 'Group Min.\n Threshold', '#keys above/equal\n threshold', '#All keys']))
	print("-"*120)
	for group in groups:
		print("-"*120)
		print( 'Analysis A not B for {}'.format(group))
		only_group_keys = groups_cKeys[group]
		if len(only_group_keys) == 0:
			print("Skipped {}...".format(group))
			print("There are no Group specific keys at specified group threshold. Change thresholds in sample details file.")
		else:
			print("-"*120)
			#print 'All keys of group above/equal to group threshold', group, '-',len(only_group_keys)
			t = []
			df_merged = pd.DataFrame(columns = ['key','#files_present','key_percent'])
			intersection_keys = []
			for group2 in groups:
				if group != group2:
					print ('Checking.. [{}] - [{}] with #group1 : {}, #group2 : {}'.format(\
						group, group2, len(only_group_keys), len(all_keys_all_groups[group2])))				
					group1_2_common = list(set(only_group_keys) & set(all_keys_all_groups[group2]))
					print ('#Common keys [{}] - [{}] : {}'.format(group, group2,len(group1_2_common)))
					if len(intersection_keys) != 0:
						intersection_keys = intersection_keys + group1_2_common
					else:
						intersection_keys = group1_2_common
					df = pd.read_csv(os.path.join(args.path, 
								args.sample_name, args.setting, 
								'keys_percents_{}_{}.csv'.format(group2, 0)),
								header = 0, index_col = 0)
					df.columns =  ['key','#files_present_'+group2, 'key_percent_'+ group2, 'total_files']
					df['key'] = df['key'].map(str)
					df['key_percent_'+ group2] = pd.to_numeric(df['key_percent_'+ group2])
					df = df[df['key'].isin(intersection_keys) ]
					# df.to_csv(os.path.join(args.path, args.sample_name, \
					# args.setting, 'intersection_keys_{}_{}.csv'.format(group, group2)))
					if len(df_merged) != 0:
						df_merged = pd.merge(df_merged, df, on = 'key', how = 'outer')
					else:
						df_merged = df
					#print(df_merged)
					only_group_keys = list(set(only_group_keys) - set(group1_2_common))
					
					print('count of remaining group A keys-', len(only_group_keys))
					t.append(len(only_group_keys))
					
				else:
					t.append(len(groups_cKeys[group]))
			t_all[group] = t
			df_merged = df_merged.fillna(0)
			df_merged.to_csv(os.path.join(args.path, args.sample_name, \
					args.setting, 'intersection_keys_without_B_filter_{}.csv'.format(group)))
			for group2 in groups:
				if group2 != group:
					df_merged = df_merged[df_merged['key_percent_'+ group2] <= float(args.b_percent)]

			df_merged.to_csv(os.path.join(args.path, args.sample_name, \
					args.setting, 'intersection_keys_{}.csv'.format(group)))
			print 'Initial only keys - ', len(only_group_keys)
			only_group_keys = only_group_keys + df_merged['key'].tolist()
			print 'after addition - ', len(only_group_keys)
			pd.DataFrame(list(only_group_keys), \
							columns = ['key'])\
							.to_excel(writer,sheet_name=get_sheet_name(group))
			summary.append((group, len(groups_cKeys[group]), len(only_group_keys)))
	print t_all
	print(summary)
	pd.DataFrame(summary, columns = [group_column_name, 'All keys in group', 'Only group keys without common'])\
						.to_excel(writer,sheet_name='summary')
	writer.close()
def get_sheet_name(name):
	if len(name)> 30:
		return name[:30]
	else:
		return name
def get_certain_percent_keys_groups_pie(percent_column_name, group_column_name):
	writer = pd.ExcelWriter(os.path.join(args.path, args.sample_name, \
					args.setting, 'common_keys_group_specific.xlsx'), \
					engine = 'xlsxwriter')
	samples_file = pd.read_csv(os.path.join(args.path, \
		args.sample_name, 'sample_details.csv'))
	groups = list(set(samples_file[group_column_name].values))
	files = glob.glob(os.path.join(args.path, args.sample_name, args.setting, '*.keys*'))
	samples_file['filename'] = os.path.join(args.path, args.sample_name, args.setting) + '//' +samples_file['protein'].astype(str).upper() + '.keys_'+args.setting
	summary = []
	groups_cKeys = {}
	all_keys_all_groups = {}
	t_all = {}
	print groups
	for group in groups:
		print group
		group_files = samples_file[samples_file[group_column_name] == group]['filename'].values
		group_threshold = list(set(samples_file[percent_column_name].values))[0]
		keys_percents = get_keys_percents(group_files, None, [group_threshold,0], group)
		groups_cKeys[group] = list(keys_percents[group_threshold])
		all_keys_all_groups[group] = list(keys_percents[0])
	for k,v in groups_cKeys.items():
		print k, '-', len(v)
	print '------------------------------------------------------------------------------------------'	
	group_intesects = {}
	for group in groups:
		only_group_keys = all_keys_all_groups[group]
		print 'All keys of group ', group, '-',len(only_group_keys)
		t = []
		intersection_keys = []
		comparing_group_count = 0
		for group2 in groups:
			if group != group2:
				comparing_group_count += 1
				print group, '-', group2
				print '------------------------------------------------------------------------------------------'
				group1_2_common = list(set(only_group_keys) & set(all_keys_all_groups[group2]))
				print group, '-', group2, 'common-',len(group1_2_common)
				df = pd.read_csv(os.path.join(args.path, args.sample_name, \
							args.setting, 'keys_percents_{}_{}.csv'.format(group2, 0)))
				df = df[df['key'].isin(group1_2_common) ]
				#df = df[df['key_percent'] <= float(args.b_percent)]
				add_keys = df['key'].values
				if len(intersection_keys) != 0:
					intersection_keys = list(set(intersection_keys) & set(add_keys))
				else:
					intersection_keys = add_keys
				group_intesects["_".join(sorted([group,group2]))] = len(group1_2_common)
				print 'Intersection of {} - {}: {}'.format(group, group2, len(group1_2_common))
				#only_group_keys = list(set(only_group_keys) - set(group1_2_common))
				
				print('count of remaining group 1 keys-', len(only_group_keys))
				t.append(len(only_group_keys))
			else:
				t.append(len(groups_cKeys[group]))
		t_all[group] = t
		print 'Initial only keys - ', len(only_group_keys)
		only_group_keys = list(only_group_keys) + list(intersection_keys)
		print 'after addition - ', len(only_group_keys)
		pd.DataFrame(list(only_group_keys), \
						columns = ['key'])\
						.to_excel(writer,sheet_name=group)
		summary.append((group, len(groups_cKeys[group]), len(only_group_keys)))
	print t_all
	pd.DataFrame(summary, columns = [group_column_name, 'All keys in group', 'Only group keys without common'])\
						.to_excel(writer,sheet_name='summary')
	writer.close()
	print group_intesects

	# Data to plot
	labels = group_intesects.keys()
	sizes = group_intesects.values()
	colors = ['gold', 'yellowgreen', 'lightcoral']
	explode = (0.1, 0, 0)  # explode 1st slice
	 
	# Plot
	plt.pie(sizes, explode=explode, labels=labels, colors=colors,
	        autopct='%1.1f%%', shadow=True, startangle=140)
	 
	plt.axis('equal')
	plt.show()

def process_files(file, req_aacds):
	print(ntpath.basename(file)[:4])
	pattern_theta_1 = ''
	pattern_dist_1 = ''
	pattern_theta_2 = ''
	pattern_dist_2 = ''
	
	samples_file = pd.read_csv(os.path.join(args.path, \
		args.sample_name, 'sample_details.csv'))
	if ntpath.basename(file)[:4] not in samples_file['protein'].values:
		return (ntpath.basename(file)[:4], None,None, None, None, None, None, None)
	if args.has_pattern:
		req_pos = re.findall(r'\d+', samples_file.set_index('protein').\
			loc[ntpath.basename(file)[:4],  'pattern'])
	else:
		req_pos = []
	records, patterns = search_aacd(open(file,'r'),req_aacds, req_pos)
	if patterns:
		# pattern_theta_1 = [pattern.split('\t')[8] for pattern in patterns]
		# pattern_dist_1 = [pattern.split('\t')[10] for pattern in patterns]
		pattern_theta_1 = patterns[0].split('\t')[8]
		pattern_dist_1 = patterns[0].split('\t')[10]
		if len(patterns) >1:
			pattern_theta_2 = patterns[1].split('\t')[8]	
			pattern_dist_2 = patterns[1].split('\t')[10]
	records.to_csv(os.path.join(args.path, args.sample_name, \
			args.setting, 'aa_search', \
	 		ntpath.basename(file)[:4] + '_'+ "_".join(req_aacds) + '.csv'))
	print (ntpath.basename(file)[:4], pd.to_numeric(records['theta']).mean(),\
	 pd.to_numeric(records['distance']).mean(), pattern_theta_1, pattern_theta_2,\
	 pattern_dist_1, pattern_dist_2, patterns)
	return (ntpath.basename(file)[:4], pd.to_numeric(records['theta']).mean(),\
	 pd.to_numeric(records['distance']).mean(), pattern_theta_1, pattern_theta_2,\
	 pattern_dist_1, pattern_dist_2, patterns)

if __name__ == '__main__':
	start = time.time()
	args = parser.parse_args()
	files = glob.glob(os.path.join(args.path,\
					 args.sample_name, args.setting, '*.triplets*'))
	column_names = ['key', 'aa0', 'pos0', 'aa1', 'pos1', 'aa2', 'pos2', \
	'classT', 'theta', 'classL', 'distance', 'x0', 'y0', 'z0', 'x1', 'y1',\
	 'z1', 'x2', 'y2', 'z2']
	
	keys_files = glob.glob(os.path.join(args.path, args.sample_name, args.setting, '*.keys*'))

	# outFolder = os.path.join(args.path, args.sample_name, args.setting, \
	# 	'common_keys')
	# if not os.path.exists(outFolder):
	# 	os.makedirs(outFolder) 
	
	# #Keys Search   
	if int(args.search_mode) == 0:
		iskeys_file = 1
		summary = []
		#Either use list of keys given in args or use common keys
		if args.use_common_keys:
			#Common Keys
			print('Calculating Common Keys')
			keys = common_keys(keys_files)
			print('Total common keys are: {}'.format(len(keys)))
			outFolder = os.path.join(args.path, args.sample_name, args.setting, \
						'common_key_search')
			pd.DataFrame(keys).to_csv(os.path.join(outFolder,'common_keys.csv'))
		else:
			if args.keys_file:
				with open(args.keys_file) as f:
					keys = f.read().splitlines()				 	
				outFolder = os.path.join(args.path, args.sample_name, args.setting, \
								'key_number_search_via_file')
			else:
				keys = [ str(key) for key in args.keys.split(',')]		
				outFolder = os.path.join(args.path, args.sample_name, args.setting, \
								'key_number_search')	
		#keys = pd.read_csv(os.path.join(outFolder,'common_keys.csv'), header = 0, index_col = 0)['key'].values
		if not os.path.exists(outFolder):
			os.makedirs(outFolder)
		for file in files:	    
			#If reading with pandas
			if iskeys_file == 1:
				df = pd.read_table(file, delimiter = '\t', names = column_names, dtype= 'unicode')
				summary.append(search_by_key(ntpath.basename(file)[:4], \
							df,\
							keys,1, outFolder))
			#line by line file read
			else:
				summary.append(search_by_key(ntpath.basename(file)[:4],\
							open(file,'r'),\
							keys, 0, outFolder))					
			
		pd.DataFrame(summary, \
			columns = ['fileName', 'all_keys', 'common_keys', 'common_keys_with_freq']) \
			.to_csv(os.path.join(outFolder,'summary_common_keys.csv'))

		calculate_low_less15_freq_from_common_keys(outFolder, args.low_freq, args.dist_less)

	# #Amino Acids Search
	if int(args.search_mode) == 2:	
		outFolder = os.path.join(args.path, args.sample_name, args.setting, \
								'aa_search')
		if not os.path.exists(outFolder):
			os.makedirs(outFolder)	
		all_list = [aas.split(',') for aas in args.aacds.upper().split(';')]
		#ploop
		#all_list = [['GLY', 'LYS','SER', 'GLY'], ['GLY', 'LYS','THR', 'GLY']]
		#lxxll
		#all_list = [['LEU', 'LEU','LEU']]
		already_completed = []
		for lst in all_list:
			print('Identifying for : {}'.format(lst))
			for req_aacds in itertools.combinations(lst, 3):
				if sorted(req_aacds) not in already_completed:
					already_completed.append(sorted(req_aacds))
					print('Started AminoAcid serach for: {}'.format("_".join(req_aacds)))
					means_theta = []
					means_distance = []
					writer = pd.ExcelWriter(os.path.join(args.path, args.sample_name, \
						args.setting, 'means_' + "_".join(req_aacds) + '.xlsx'), \
						engine = 'xlsxwriter')
					
					means_theta_distance  = Parallel(n_jobs=cpu_count() - 1, verbose=10, \
						backend="multiprocessing", batch_size="auto")(\
						delayed(process_files)(file, req_aacds) for file in files)
					pd.DataFrame(means_theta_distance, \
						columns = ['file', 'all_theta_per_protein', \
								'all_dist_per_protein', 'motif_theta_1','motif_theta_2', \
								'motif_dist_1', 'motif_dist_2','motif'])\
						.to_excel(writer,sheet_name='means_theta_distance')
					writer.close()

	#Sequence Identification
	if int(args.search_mode) == 3:
		writer_pattern = pd.ExcelWriter(os.path.join(args.path, args.sample_name, \
					args.setting, 'patterns_{}{}'.format(args.identifying_pattern, '.xlsx')), \
					engine = 'xlsxwriter')
		#ploop
		if args.identifying_pattern == 'ploop':
			patterns = [['GLY','GLY','LYS','SER'], ['GLY','GLY','LYS','THR']]
			no_of_xs = 4
		if (args.identifying_pattern == 'leu1') or (args.identifying_pattern == 'leu2'):
			patterns = [['LEU', 'LEU','LEU']]
			no_of_xs = 2
		if args.identifying_pattern == 'gen':
			patterns = [args.gen_pattern_aacds.upper().split(',')]
			no_of_xs = int(args.no_of_xs)
		for pattern in patterns:
			matched_list = []
			# sequences  = Parallel(n_jobs=cpu_count() - 1, verbose=10, \
			# 			backend="multiprocessing", batch_size="auto")(\
			# 			delayed(identify_sequence)(file, column_names, \
			#  	[r'GLY_\d+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+GLY_\d+LYS_\d+SER_\d+', \
			#  	r'GLY_\d+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+GLY_\d+LYS_\d+THR_\d+']) for file in files)
			# identify_sequence(files,column_names, \
			# 	[r'GLY_\d+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+GLY_\d+LYS_\d+SER_\d+', \
			# 	r'GLY_\d+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+GLY_\d+LYS_\d+THR_\d+'])

			#Uncomment this for parallel execution
			# matched_list = Parallel(n_jobs=cpu_count() - 1, verbose=10, \
			# 			backend="multiprocessing", batch_size="auto")(\
			# 			delayed(identify_pattern_by_normal)(file,pattern, no_of_xs, args.identifying_pattern)\
			# 			for file in files)
			for file in files:
				matched_list.append(identify_pattern_by_normal(file,pattern, no_of_xs, args.identifying_pattern))
			# for file in files:
			# 	matched_list.append(identify_pattern_by_normal(file,pattern, no_of_xs))
			pd.DataFrame(matched_list, columns = ['file','pattern_matched'])\
				.to_excel(writer_pattern, sheet_name = "_".join(pattern))

	if int(args.search_mode) == 4:
		get_common_keys_groups(args.get_common_keys_groups)
	if int(args.search_mode) == 5:
		outFolder = os.path.join(args.path, args.sample_name, args.setting)
		calculate_low_less15_freq_from_common_keys(outFolder, args.low_freq, args.dist_less)
	if int(args.search_mode) == 6:
		print("Entered mode 6")
		files = glob.glob(os.path.join(args.path,\
					 args.sample_name, args.setting, '*.keys_*'))
		get_certain_percent_keys_groups(args.percent_column_name, args.group_column_name)
		#get_certain_percent_keys_groups_pie()

	end = time.time()
	print("Task Completed. Total time taken: {} mins".format((end - start)/60))	

	

	
        


