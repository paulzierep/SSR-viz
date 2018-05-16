import datetime
from matplotlib.backends.backend_pdf import PdfPages
from .lib.mpl_strong_colors import cnames

from scipy import stats
import math
import pandas as pd
import numpy as np
import string
import itertools
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
#for sub_matrix if the matrix is not provided
from Bio.Data import IUPACData
import Bio.SubsMat.MatrixInfo as MI

from . import sequence_translations as st
#import random  
#sys.path.append('../')
from .basic import basic
#import basic
import os

################################
# basic functions for the plots
################################

def label_from_ind_plot(ind):
	label = ' vs '.join([ind[0], ind[1]])
	if ind[2]:
		sign = ' (w)'
	else:
		sign = ''

	label += sign
	return(label)

def label_from_ind_heatmap(ind):
	label = ' vs '.join([ind[0], ind[1]])
	return(label)

def multi_pdf(plotlist, path = 'temp.pdf', info = ('PDF', '-', '-', '-')):

	if '.pdf' not in path:
		path += '.pdf'

	with PdfPages(path) as pdf:


		index = 1
		for fig in plotlist:
			#print(fig)
			#plt.title(str(index))
			pdf.savefig(fig)  # saves the current figure into a pdf page
			#plt.close()
			index += 1

		# We can also set the file's metadata via the PdfPages object:
		d = pdf.infodict()
		d['Title'] = info[0]
		d['Author'] = info[1]
		d['Subject'] = info[2]
		d['Keywords'] = info[3]
		d['CreationDate'] = datetime.datetime.today()


def cm2inch(value):
	return(value/2.54)

# input - df: a Dataframe, chunkSize: the chunk size
# output - a list of DataFrame
# purpose - splits the DataFrame into smaller of max size chunkSize (last is smaller)
def splitDataFrameIntoSmaller_col(df, chunkSize = 200): 
	
	# print(df)
	# exit()

	listOfDf = list()
	numberChunks = len(df.columns) // chunkSize + 1
	#print(numberChunks)

	for i in range(numberChunks):
		begin_idx = i*chunkSize
		if i != 0:
			begin_idx += 1
		end_idx = (i+1)*chunkSize

		#print(df.loc[:,begin_idx:end_idx])
		#print(df.loc[:,begin_idx:end_idx])
		df_slice = df.loc[:,begin_idx:end_idx]
		if not df_slice.empty:
			listOfDf.append(df_slice)
	
	return(listOfDf)

def chunck_df_col(df, chunksize = 200): 

	df.columns = df.columns.map(int)

	df_list = []
	# begin_index = 
	# print(chunksize)
	# exit()
	for chunk in np.arange(df.columns[0], df.columns[-1], chunksize):
		#print(chunk)
		begin = chunk
		end = chunk + chunksize - 1
		if end > df.columns[-1]:
			end = df.columns[-1]

		new_df = df.loc[:,begin:end]

		df_list.append(new_df)

	return(df_list)

def export2jalview(df, annot = '', path = None):

	colors = itertools.cycle(iter(list(cnames.items()))) #automatic colors

	init_text = 'JALVIEW_ANNOTATION \n'

	df = df.round(2)
	df[df < 0 ] = 0
	for ind in df.index:
		matrix = list(df.loc[ind].as_matrix())
		matrix_2 = [str(x) + ',' + str(int(round(x*10,0))) for x in matrix]
		matrix_text = '|'.join(matrix_2)

		jal_line = 'BAR_GRAPH\t{0}\tSSP {1} \t{2}\n'.format(label_from_ind_plot(ind), annot,matrix_text)
		init_text += jal_line
		jal_line2 = 'COLOUR\t{0}\t{1}\n'.format(label_from_ind_plot(ind), colors.next()[1])
		init_text += jal_line2

	with open(path, 'w') as jal_file:
		jal_file.write(init_text)


#######################################
# cons diff algos
#######################################


def get_sub_matrix(name = 'blosum100', gap_importance = 1):

	#get an alphabet
	p_letters_gap = IUPACData.extended_protein_letters + '-'

	if name == 'basic':
		#create basic matrix
		sub_matrix = pd.DataFrame(np.identity(len(p_letters_gap)),\
					index=list(p_letters_gap), \
					columns=list(p_letters_gap))

	else:
		#get an matrix
		p_letters_gap = IUPACData.extended_protein_letters + '-'
		sub_matrix = pd.DataFrame(np.zeros((len(p_letters_gap),len(p_letters_gap))),
						index=list(p_letters_gap), \
						columns=list(p_letters_gap))

		if name in MI.available_matrices:
			matrix = getattr(MI,name)
		else:
			print(('Matrix {0} does not exsist, chosse one of those {1}'.format(name, MI.available_matrices)))

		for AS, score in list(matrix.items()):
			sub_matrix.loc[AS[0],AS[1]] = score
			sub_matrix.loc[AS[1],AS[0]] = score

	#normalize the matrix, sot that it works with the algorithm
	sub_matrix = sub_matrix.apply(lambda row: row + -min(row), axis = 0)
	sub_matrix = sub_matrix.apply(lambda row: 1 - (row / max(row)), axis = 0)

	#special ASs are not used at the moment
	sub_matrix.loc[['U','O','J'],:] = 0
	sub_matrix.loc[:,['U','O','J']] = 0

	#the importance of the gap can be set individually from the matrix
	sub_matrix.loc['-',:] = gap_importance
	sub_matrix.loc[:,'-'] = gap_importance
	sub_matrix.loc['-','-'] = 0

	#some matrices do not support all letters
	sub_matrix.fillna(value = 0, inplace = True)

	return(sub_matrix)

def conservation_scores(df1, df2, sub_matrix = None, cons_type = 'entropy'):
	'''
	Calculates the conservation scores between df1 and df2 
	(should be two arrays of identical dimensions), 
	if df = n*m --> sub_matrix should be n*n, 
	returns a 1d array with the cons scores
	'''

	if not isinstance(sub_matrix, pd.DataFrame):
		print('Sub matrix must be a Dataframe !')
		exit()

	# print(df1.loc[2].to_dict())
	# exit()

	df1 = df1.reindex(sub_matrix.columns, axis = 1)
	df2 = df2.reindex(sub_matrix.columns, axis = 1)

	#the calculation of the croneberg product only works with numpy so far for me
	#but the indices must correspond to the correct column

	sub_matrix = np.array(sub_matrix)

	df1 = np.array(df1)                                                         #convert to np array
	df2 = np.array(df2)

	#########################
	# conservedness 
	#########################

	########################
	# implement entropy !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	########################

	cons1 = np.max(df1, axis = 1)
	# print(cons1)

	if cons_type == 'entropy':
		#normalized inverted shanon entropy
		entro_func = lambda x: 1 - stats.entropy(x, qk = None, base = float(2)) / 4.32

		cons1 = np.apply_along_axis(entro_func, 1, df1)
		cons2 = np.apply_along_axis(entro_func, 1, df2)

	else:

		cons1 = np.max(df1, axis = 1)                                               #cons score can be defined as the
		cons2 = np.max(df2, axis = 1)

	cons_score = (cons1 + cons2) / 2                                            #both sequences are considered, normalized to 1




	#########################
	# entropy
	#########################


	sub_probability = (df1[...,None]*df2[:,None,:])                             #multiply row by row, see https://stackoverflow.com/questions/35162197/numpy-matrix-multiplication-of-2d-matrix-to-give-3d-matrix
																				#--> substitution probaility matrix

	sub_diff_matrix = sub_probability*sub_matrix                                # only different AS are interesting 


	sub_diff_score = np.sum(np.sum(sub_diff_matrix, axis=1), axis=1)            #the actual differnet score of each row

	cons_diff_score = cons_score * sub_diff_score                               #final score considering cons. in a pssm
																				#and differences between them 


	return(cons_diff_score)


def conservation_scores_2(df1, df2, sub_matrix = None):
	'''
	Calculates the conservation scores between df1 and df2 
	(should be two arrays of identical dimensions), 
	if df = n*m --> sub_matrix should be n*n, 
	returns a 1d array with the cons scores
	'''

	if not isinstance(sub_matrix, pd.DataFrame):
		print('Sub matrix must be a Dataframe !')
		exit()
	# 	#length = len(df1[])
	# 	print('There is something todo here')
	# 	exit()
	# 	sub_matrix = pd.DataFrame(np.identity(length),\
	# 					index=list(string.ascii_uppercase[:length]), \
	# 					columns=list(string.ascii_uppercase[:length]))      #identity matrix --> substitutuion m. like Blossom...       


	# print(sub_matrix)
	# exit()



	# print(sub_matrix)
	# print(df1.reindex_axis(sub_matrix.columns, axis = 1))


	df1#.reindex_axis(sub_matrix.columns, axis = 1)
	df2#.reindex_axis(sub_matrix.columns, axis = 1)

	print(df1)
	exit()
	#print()

	#sub_matrix = np.array(sub_matrix)                                         #inverse matrix where 0 = 1, 1 = 0



	#print(sub_matrix)

	df1 = np.array(df1)                                                         #convert to np array
	df2 = np.array(df2)

	#########################
	# conservedness 
	#########################

	########################
	# implement entropy !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	########################

	cons1 = np.max(df1, axis = 1)                                               #cons score can be defined as the 
																				#max of the row in the pssm
	cons2 = np.max(df2, axis = 1)
	cons_score = (cons1 + cons2) / 2                                            #both sequences are considered, normalized to 1


	#########################
	# entropy
	#########################


	sub_probability = (df1[...,None]*df2[:,None,:])                             #multiply row by row, see https://stackoverflow.com/questions/35162197/numpy-matrix-multiplication-of-2d-matrix-to-give-3d-matrix
																				#--> substitution probaility matrix

	sub_diff_matrix = sub_probability*sub_matrix                                # only different AS are interesting 


	sub_diff_score = np.sum(np.sum(sub_diff_matrix, axis=1), axis=1)            #the actual conservation score of each row

	cons_diff_score = cons_score * sub_diff_score                               #final score considering cons. in a pssm
																				#and differences between them 


	return(cons_diff_score)

'''
automatic creation of test data, for conservation_scores function test
'''

# size = 100                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
# length = 26 #max 26

# sub_matrix = pd.DataFrame(np.identity(length),\
#                       index=list(string.ascii_uppercase[:length]), \
#                       columns=list(string.ascii_uppercase[:length]))      #identity matrix --> substitutuion m. like Blossom...       

# df1 = pd.DataFrame(np.random.dirichlet(np.ones(length),size=size),\
#                       columns=list(string.ascii_uppercase[:length]))      #random df representing a pssm 
#                                                                           # dirichlet(np.ones(26) -> (random numbers that sum up to 1
#                                                                           # string.ascii_uppercase -> all letters in alphabet (ascii)
# df2 = pd.DataFrame(np.random.dirichlet(np.ones(length),size=size), \
#                       columns=list(string.ascii_uppercase[:length]))      #random df representing another pssm

# #examples

###super diff
# df1.loc[0] = 0    #set all to 0
# df1['A'][0] = 1   #set one to 1

# df2.loc[0] = 0    #set all to 0
# df2['B'][0] = 1   #set one to 1

##super diff distributed
# df1.loc[0] = 0        #set all to 0
# df1['A'][0] = 0.5     #set one to 1
# df1['B'][0] = 0.5     #set one to 1

# df2.loc[0] = 0        #set all to 0
# df2['C'][0] = 0.5     #set one to 1
# df2['D'][0] = 0.5     #set one to 1

##super diff distributed
# df1.loc[0] = 0            #set all to 0
# df1['A'][0] = 0.3         #set one to 1
# df1['B'][0] = 0.3         #set one to 1

# df2.loc[0] = 0            #set all to 0
# df2['C'][0] = 0.3         #set one to 1
# df2['D'][0] = 0.3         #set one to 1


# print(df1)
# print(df2)
# print(conservation_scores(df1, df2, sub_matrix))

# exit()

###########################
# old version
###########################

# def all_vs_all_pssm_cons(pssm_dict, sub_matrix = pd.DataFrame()):
#   '''
#   calls the conservation_scores function on all 
#   class combinations, retuns a pandas df, with all the results, 
#   positions as columns, class as rows (doublicated) -> 
#   allows easy use of the df as matrix 
#   '''
#   row = []

#   #create all combinations (no doublicates) using itertools
#   discri_keys = pssm_dict.keys()
#   combinations = list(itertools.combinations(discri_keys, 2)) #Example: ('3', '2'), ('3', '5'), ('3', '4')
	
#   #create an empty array of the expected result file -> combis * positions
#   #combi times tow -> seems to be easier to create a df with doublicated rows, then double indices 
#   num_combinations = len(combinations)
#   num_positions = pssm_dict[pssm_dict.keys()[0]].shape[0]
#   #result_array = np.zeros((num_combinations, num_positions))
#   result_array = np.zeros((2*num_combinations, num_positions))

#   #call the conservation_scores for each pssm and write into result table
#   index = 0
#   for combi in combinations:
#       row.append(combi[0])
#       row.append(combi[1])
#       single_cons_score = conservation_scores(pssm_dict[combi[0]], pssm_dict[combi[1]], sub_matrix)
#       result_array[index] = single_cons_score
#       result_array[index + 1] = single_cons_score
#       index += 2
#       #index += 1

#   df = pd.DataFrame(result_array, index = row)
#   return(df)

def all_vs_all_pssm_cons(pssm_dict, sub_matrix = pd.DataFrame()):
	'''
	calls the conservation_scores function on all 
	class combinations, retuns a pandas df, with all the results, 
	positions as columns, class as rows (doublicated) -> 
	allows easy use of the df as matrix 
	'''

	if sub_matrix.empty:
		p_letters_gap = IUPACData.extended_protein_letters + '-'
		sub_matrix = pd.DataFrame(np.identity(len(p_letters_gap)),\
						index=list(p_letters_gap), \
						columns=list(p_letters_gap))        #identity matrix --> substitutuion m. like Blossom...       


	row = []

	#create all combinations (no doublicates) using itertools
	discri_keys = list(pssm_dict.keys())
	combinations = list(itertools.combinations(discri_keys, 2)) #Example: ('3', '2'), ('3', '5'), ('3', '4')
	
	#create an empty array of the expected result file -> combis * positions
	#combi times tow -> seems to be easier to create a df with doublicated rows, then double indices 
	num_combinations = len(combinations)
	num_positions = pssm_dict[list(pssm_dict.keys())[0]].shape[0]
	result_array = np.zeros((num_combinations, num_positions))
	#result_array = np.zeros((2*num_combinations, num_positions))

	#call the conservation_scores for each pssm and write into result table
	index = 0
	c1 = []
	c2 = []
	for combi in combinations:
		c1.append(combi[0])
		c2.append(combi[1])
		#row.append(combi[1])
		single_cons_score = conservation_scores(pssm_dict[combi[0]], pssm_dict[combi[1]], sub_matrix)
		result_array[index] = single_cons_score
		#result_array[index + 1] = single_cons_score
		#index += 2
		index += 1

	df = pd.DataFrame(result_array)
	df['Class_A'] = c1
	df['Class_B'] = c2

	return(df)

def add_all_vs_all(df):

	'''
	computes the average cons_score for the all_vs_all matrix
	'''

	#total all_vs_all

	df_ref = df.drop(['Class_A', 'Class_B'], axis = 1)
	df_total_sum = df_ref.mean(axis=0)
	df_total_sum = pd.DataFrame([df_total_sum], columns = df.columns)
	df_total_sum['Class_A'] = 'all'
	df_total_sum['Class_B'] = 'all'


	#one vs all 
	one_vs_all = []
	for clas in set(list(df['Class_A']) + list(df['Class_B'])):
		df_ref = df.drop(df[(df['Class_A'] != clas) & (df['Class_B'] != clas)].index)
		df_sum = df_ref.mean(axis=0)
		df_sum = pd.DataFrame([df_sum], columns = df.columns)
		df_sum['Class_A'] = clas
		df_sum['Class_B'] = 'all'
		one_vs_all.append(df_sum)

	final_df = pd.concat([df,df_total_sum] + one_vs_all, ignore_index=True)
	final_df.set_index(['Class_A', 'Class_B'], inplace = True)

	return(final_df)


def df2pssm_visual(df, path=None,
				 **kwargs):


	#most alignemnt tools start with 1 instead of 0, so this can be changed
	if 'alignment_index' in kwargs.keys():
		if kwargs['alignment_index']:
			df.columns = df.columns + 1

	#print(df.loc[('all', 'all'),305:320])

	df.sort_index(inplace = True)

	###########################
	# Arguments for the heatmap
	###########################

	#generate a choice of all possible classes, without all
	optional_classes = list(set(list(df.index.levels[1]) + list(df.index.levels[0])))
	optional_classes.remove('all')

	all_df_list = []

	#get the all vs all dataframes
	if basic.check_dict_key_not_none(kwargs, 'hm_all'):

		if kwargs['hm_all'] == []: #if not specified all are taken
			kwargs['hm_all'] = optional_classes

		df_na = df[(df.index.get_level_values('Class_A') != 'all')]
		df_na = df_na[(df_na.index.get_level_values('Class_B') != 'all')]

		class2drop = [x for x in optional_classes if x not in kwargs['hm_all']] #drop all classes, which are not present 

		for clas in class2drop:
			row2drop = df_na.loc[(df_na.index.get_level_values('Class_A') == clas) | (df_na.index.get_level_values('Class_B') == clas)].index
			df_na = df_na.drop(row2drop)

		all_df_list.append(df_na)

	#get the one vs all dataframes
	if basic.check_dict_key_not_none(kwargs, 'hm_ova'):
		if kwargs['hm_ova'] == []:
			kwargs['hm_ova'] = optional_classes

		df_na = df[(df.index.get_level_values('Class_A') != 'all')]
		df_na = df_na[(df_na.index.get_level_values('Class_B') == 'all')]

		for clas in kwargs['hm_ova']:

			hm_all_df = df_na.loc[(df_na.index.get_level_values('Class_A') == clas)]
			all_df_list.append(hm_all_df)

		# #get the all vs all dataframe
		# print(kwargs['hm_ava'])
	if not kwargs['no_hm_ava']:
		df_na = df[(df.index.get_level_values('Class_A') == 'all')]
		all_df_list.append(df_na)

	if not all_df_list:
		print('''There are no values given for the heatmap, 
			please add at least the all-vs-all representation by unchecking the 
			hm_ava tick
			''')
		exit()

	else:
		df_heatmap = pd.concat(all_df_list)

	###########################
	# Arguments for the plot
	###########################


	all_df_list = []

	#get the all vs all dataframes
	if basic.check_dict_key_not_none(kwargs, 'pl_all'):
		if kwargs['pl_all'] == []: #if not specified all are taken
			kwargs['pl_all'] = optional_classes

		df_na = df[(df.index.get_level_values('Class_A') != 'all')]
		df_na = df_na[(df_na.index.get_level_values('Class_B') != 'all')]

		class2drop = [x for x in optional_classes if x not in kwargs['pl_all']]

		for clas in class2drop:
			row2drop = df_na.loc[(df_na.index.get_level_values('Class_A') == clas) | (df_na.index.get_level_values('Class_B') == clas)].index
			df_na = df_na.drop(row2drop)

		all_df_list.append(df_na)

	#get the one vs all dataframes

	if basic.check_dict_key_not_none(kwargs, 'pl_ova'):
		if kwargs['pl_ova'] == []:
			kwargs['pl_ova'] = optional_classes

		df_na = df[(df.index.get_level_values('Class_A') != 'all')]
		df_na = df_na[(df_na.index.get_level_values('Class_B') == 'all')]

		for clas in kwargs['pl_ova']:
			#df_na = df[(df.index.get_level_values('Class_B') == 'all')]
			hm_all_df = df_na.loc[(df_na.index.get_level_values('Class_A') == clas)]
			all_df_list.append(hm_all_df)

	#get the all vs all dataframe
	if not kwargs['no_pl_ava']:
		df_na = df[(df.index.get_level_values('Class_A') == 'all')]
		all_df_list.append(df_na)

	if not all_df_list:
		print('''There are no values given for the plot, 
			please add at least the all-vs-all representation by unchecking the 
			pl_ava tick
			''')
		exit()

	else:
		df_plot = pd.concat(all_df_list)

	###########################
	# Special arguments for the plot
	###########################

	if basic.check_dict_key_not_none(kwargs, 'window'):
		window =  int(kwargs['window'])

		if basic.check_dict_key_not_none(kwargs, 'window_type'):
			window_type =  kwargs['window_type']
		else:
			window_type = 'mean'

		if window_type == 'mean':

			df_roll = df_plot.rolling(window, \
				center= True, min_periods=1, axis = 1).mean()          #apply window function to exclude only long period of gaps

		if window_type == 'max':
			df_roll = df_plot.rolling(window, \
				center= True, min_periods=1, axis = 1).max()          #apply window function to exclude only long period of gaps

		if window_type == 'min':
			df_roll = df_plot.rolling(window, \
				center= True, min_periods=1, axis = 1).min()          #apply window function to exclude only long period of gaps

		if window_type == 'std':
			df_roll = df_plot.rolling(window, \
				center= True, min_periods=1, axis = 1).std()          #apply window function to exclude only long period of gaps

		df_roll['Window'] = True
		df_roll.set_index('Window', append=True, inplace=True)

		df_plot['Window'] = False
		df_plot.set_index('Window', append=True, inplace=True)
		df_plot = pd.concat([df_plot, df_roll])

	else:

		df_plot['Window'] = False
		df_plot.set_index('Window', append=True, inplace=True)


	def drop_lowest(row, num):
		if row.name[2]:
			return(row)
		
		drop_panalty = min(row.nlargest(num)) #get a panaly which allows to drop the n below that panelty
		columns_to_drop = row < drop_panalty  #get df where condition applies
		row[columns_to_drop] = -1         #set to nan
		return(row)

	def drop_by_panalty(row, panalty):
		if not row.name[2]:
			return(row)

		columns_to_drop = row < panalty  #get df where condition applies
		row[columns_to_drop] = -1         #set to nan
		return(row)

	if 'get_best' in kwargs:
		df_plot = df_plot.apply(lambda row: drop_lowest(row, kwargs['get_best']), axis =1)
		#df_plot.apply(lambda row: drop_lowest(row, kwargs['get_best']), axis =1)
		# print(df_plot)

	if 'drop_panalty' in kwargs:
		df_plot = df_plot.apply(lambda row: drop_by_panalty(row, kwargs['drop_panalty']), axis =1)

	###########################
	# Plot options
	###########################

	if basic.check_dict_key_not_none(kwargs, 'figsize'):
		figsize = (float(kwargs['figsize'][0]), float(kwargs['figsize'][1]))
	else:
		figsize = (29.7, 21.0)

	#print(figsize)
	# print(kwargs['chunksize'])

	if  basic.check_dict_key_not_none(kwargs, 'chunksize'):
		
		if kwargs['chunksize'] == 'total':
			chunksize = len(df.columns)
		else:
			chunksize = int(kwargs['chunksize'])
	else:
		chunksize = 100

	if basic.check_dict_key_not_none(kwargs, 'fontsize'):
		fontsize = kwargs['fontsize']
	else:
		fontsize = 8

	if basic.check_dict_key_not_none(kwargs, 'w_ratio'):
		w_ratio = (float(kwargs['w_ratio'][0]), float(kwargs['w_ratio'][1]))
	else:
		w_ratio = (1, 3)

	if basic.check_dict_key_not_none(kwargs, 'tick_ratio'):
		tick_ratio = int(kwargs['tick_ratio'])
	else:
		tick_ratio = 5

	######################
	# Create plots for Jalview
	######################

	if basic.check_dict_key_true(kwargs, 'jv_plot'):
		jv_path = path + '_jv_plot.txt'
		export2jalview(df_plot, annot = 'plot', path = jv_path)

	# if basic.check_dict_key_true(kwargs, 'jv_heatmap'):
	# 	jv_path = path + '_jv_heatmap.txt'
	# 	export2jalview(df_heatmap, annot = 'heatmap', path = jv_path)

	#################################
	# Plotting
	#################################

	listofdf_plot = chunck_df_col(df_plot, chunksize = chunksize)
	listofdf_heatmap = chunck_df_col(df_heatmap, chunksize = chunksize)

	#get one fig per df
	figure_list = []
	index = 0

	for df_plot, df_heatmap in zip(listofdf_plot, listofdf_heatmap):

		df_plot.sort_index(axis = 0, level = 0, inplace = True)
		df_heatmap.sort_index(axis = 0, level = 0, ascending=False ,inplace = True)

		#print(df_heatmap)

		#get a fig object, an ax for the plot, another for the heatmap
		if basic.check_dict_key_true(kwargs, 'top_label'):
			fig, (ax_label, ax, ax2) = plt.subplots(3,1, gridspec_kw = {'height_ratios':[0.25,w_ratio[0],w_ratio[1]]}, figsize=(cm2inch(figsize[0]), cm2inch(figsize[1])))
		else:
			fig, (ax, ax2) = plt.subplots(2,1, gridspec_kw = {'height_ratios':[w_ratio[0],w_ratio[1]]}, figsize=(cm2inch(figsize[0]), cm2inch(figsize[1])))

		##################
		# top plot
		##################

		colors = itertools.cycle(iter(list(cnames.items()))) #automatic colors

		# if basic.check_dict_key_true(kwargs, 'plot_all'): #all plot option
		#print(df_plot.sort_index(axis = 0, level=0))
		# print(df_plot)
		# exit()
		# exit()
		plot_labels = []
		color = next(colors)[1]
		for ind in df_plot.index:

			label = label_from_ind_plot(ind)

			if ind[2]:
				ax.plot(df_plot.columns, df_plot.loc[ind], '-' ,markersize=8,label = ind, color = color)
			else:
				# print(df_plot.columns)
				# print(df_plot.loc[ind])
				color = next(colors)[1]

				#print(ind)
				if ind[0] == 'all' and ind[1] == 'all':
					marker = 'D'
				elif ind[0] != 'all' and ind[1] == 'all':
					marker = 'v'
				else:
					marker = 'o'


				ax.plot(df_plot.columns, df_plot.loc[ind],marker,markersize=8,label = ind, color = color)

			plot_labels.append(label)


		# ax.plot(df_plot.columns, df_plot.loc['Average'], 'o', markersize=8, label = 'Average', color = 'blue')

		handles, labels = ax.get_legend_handles_labels()

		if basic.check_dict_key_true(kwargs, 'top_label'):
			ax_label.legend(handles, plot_labels, loc='upper center', ncol = 7, fontsize = fontsize)
			ax_label.axis('off')


		x_ticks_major = np.arange(min(df_plot.columns), max(df_plot.columns)+ 1, tick_ratio)
		x_ticks_minor = np.arange(min(df_plot.columns), max(df_plot.columns)+ 1, 1)

		ax.set_xticks(x_ticks_major)
		ax.set_xticks(x_ticks_minor, minor=True)
		ax.set_xticklabels(x_ticks_major, fontsize = fontsize)
		ax.set_ylim(0 , 1)
		ax.set_yticklabels([0,0.2,0.4,0.6,0.8,1], fontsize = fontsize)

		min_x = min(df_plot.columns) - 1
		max_x = (index + 1) * chunksize + 1


		ax.set_xlim(min_x,max_x)
		ax.set_ylabel('Score', fontsize = fontsize)

		#####################
		# heatmap below
		#####################

		#df.columns = df.columns + 0.5

		ax2.pcolor(df_heatmap, cmap='YlOrRd', vmin = 0, vmax = 1)

		ticks_range = np.arange(0.5, len(df_heatmap.index), 1)

		if basic.check_dict_key_true(kwargs, 'drop_class_label'):
			ax2.set_yticks([np.mean(ticks_range)])
			ax2.set_yticklabels(['-'], fontsize = fontsize)

		else:
			ax2.set_yticks(ticks_range)
			labels = [label_from_ind_heatmap(ind) for ind in df_heatmap.index]
			ax2.set_yticklabels(labels, fontsize = fontsize)

		# print(ax2.axis())
		# print(df.columns)

		min_x = - 0.5
		max_x = chunksize + 0.5
		ax2.set_xlim(min_x,max_x)

		x_ticks_major_ax2 = x_ticks_major - (index*chunksize) - 0.5
		ax2.set_xticks(x_ticks_major_ax2)
		ax2.set_xticklabels(x_ticks_major, fontsize = fontsize)

		ax2.set_xlabel('Position in alignment', fontsize = fontsize)
		ax2.set_ylabel('Class labels', fontsize = fontsize)

		fig.subplots_adjust(hspace=0.1)

		figure_list.append(fig)

		print(('fig {0} created'.format(index)))

		index += 1


	fig, (ax, ax2) = plt.subplots(2,1, gridspec_kw = {'height_ratios':[1, 4]}, figsize=(cm2inch(figsize[0]), cm2inch(figsize[1])))

	fig.subplots_adjust(top = 0.6, hspace = 0.3)

	ax2.legend(handles, plot_labels, loc='upper center', ncol = 7, fontsize = fontsize)
	ax2.axis('off')

	#create a options text for the top plot
	kwargs2scip = ['keep_data_folder','no_rule','no_alignment','visual','command_line']
	kwargs_txt = ''
	i = 1
	for key, item in list(kwargs.items()):
		if item and key not in kwargs2scip:
			kwargs_txt += key + ' : ' + str(item) + '    '
			if i % 4 == 0:
				kwargs_txt += '\n'
			i += 1

	# if basic.check_dict_key_not_none(kwargs, 'command_line'):
	# 	command_line = kwargs['command_line']
	# else:
	# 	command_line = 'Not supplied'

	#fig.suptitle('Options used for this plot:\n{0}\n Command line input: {1}'.format(kwargs_txt, command_line),
	fig.suptitle('Options used for this plot:\n{0}'.format(kwargs_txt),
	fontsize = fontsize, fontweight='bold')

	cb1 = mpl.colorbar.ColorbarBase(ax, cmap='YlOrRd', orientation='horizontal')
	cb1.ax.tick_params(labelsize=fontsize)
	cb1.set_label('Colorbar range', fontsize = fontsize)

	figure_list = [fig] + figure_list

	multi_pdf(figure_list, path = path)



def slice_function(name, class_name):
	'''
	Return True if the class_name is part of the name based on the pattern (bla vs blu)
	'''
	name = name.split('vs')
	name = [item.strip(' ') for item in name]

	if class_name in name:
		return(True)
	else:
		return(False)


def slice_class_of_df(df, class_name):
	'''
	Slices only one class from the df
	'''
	inds = [ind for ind in df.index if slice_function(ind, class_name)]
	df = df.loc[inds]
	return(df)
	



#####################################
# Scoring of the pssm algo
#####################################


def pssm_scoring(scoring_seq, pssm_dict):

	scoring_seq_pssm = st.seq_list2pssm([scoring_seq])                              #create a pssm of the scoring sequence (only 1 or 0)

	scoring_dict = {}
	for clas, pssm in list(pssm_dict.items()):
		scored_pssm = pssm * scoring_seq_pssm                                       #multiply df to obtain scored pssm
		scored_pssm = scored_pssm.sum(axis = 1).to_frame()                          #get sum of all rows -> actual pssm score
		scored_pssm = scored_pssm.T                                                 #better for next function
		scoring_dict[clas] = scored_pssm

	return(scoring_dict)


def weighed_pssm_scoring(scoring_dict, weigth_matrix, no_weigth_matrix = True, position_threshold = None):

	result_dict = {}

	# plain pssm based scoring (sum of all the scores)
	if no_weigth_matrix:
		for clas, scores in list(scoring_dict.items()):
			result_dict[clas] = float(scores.sum(axis = 1))

		return(result_dict)

	if not position_threshold:
		print('No position threshold -> entire alignment is used')
		position_threshold = len(weigth_matrix.columns)

	if position_threshold > len(weigth_matrix.columns):
		print('The position threshold must be smaller then the alignment -> entire alignment is used')
		position_threshold = len(weigth_matrix.columns)

	# the modified weight matrix

	####################
	#todo add more modification possibilities
	####################

	mod_weigth_matrix = weigth_matrix.apply(lambda x: x.sort_values(ascending=False).head(position_threshold), axis=1) # <---works

	#create empty result table, to store the info
	result_table = pd.DataFrame(index = mod_weigth_matrix.index, columns = list(scoring_dict.keys()))

	#get the labes from the result table
	lables = mod_weigth_matrix.index

	# fill the result table, a bit messy due to 'vs' handling, works anyway

	for clas, scores in list(scoring_dict.items()):

		# print clas
		# print scores
		# print lables

		for lable in lables:
			s_lable = lable.split(' vs ')
			if clas == s_lable[0] or clas == s_lable[1]:

				weighed_scoring = scores * mod_weigth_matrix.loc[lable]
				weighed_scoring = weighed_scoring.sum(axis = 1)

				result_table.loc[lable, clas] = float(weighed_scoring)

	result_table = add_average(result_table)

	####################
	#todo figure of result table
	####################

	#result table to result result_dict

	result_dict = result_table.loc['Average'].to_dict()

	return(result_dict)


###################################
# unit test
###################################

'''
automatic creation of test data for the all-vs-all pssm conservation_scores function
'''

# import sys
# import os
# sys.path.append('../')

# import sequence_translations as st

# from Bio.Data import IUPACData

# import domain_extractor

# # ###############
# # # create the sub_matrix
# # ###############

# p_letters_gap = IUPACData.extended_protein_letters + '-'

# sub_matrix = pd.DataFrame(np.identity(len(p_letters_gap)),\
#                       index=list(p_letters_gap), \
#                       columns=list(p_letters_gap))                                #identity matrix --> substitutuion m. like Blossom...       

# ###############
# # get test seq and discri dict
# ###############

# length = 260

# discri_dict = st.get_random_discri_dict(num_classes = 10, num_seqs = 10, length = length, seed_diff = 0.1, class_diff = 0.2, random_seed = True)
# #test_seq = st.get_random_seqs(num = 1, length = length, difference = 0.6, initial_seq = None, names = None)

# # # # # ###############
# # # # # # compute weigth_matrix
# # # # # ###############

# pssm_dict = st.discri_dict2pssm_dict(discri_dict) 

# '''
# ['benner6', 'benner22', 'benner74', 'blosum100',
# 'blosum30', 'blosum35', 'blosum40', 'blosum45',
# 'blosum50', 'blosum55', 'blosum60', 'blosum62',
# 'blosum65', 'blosum70', 'blosum75', 'blosum80',
# 'blosum85', 'blosum90', 'blosum95', 'feng', 'fitch',
# 'genetic', 'gonnet', 'grant', 'ident', 'johnson', 'levin',
# 'mclach', 'miyata', 'nwsgappep', 'pam120', 'pam180', 'pam250',
# 'pam30', 'pam300', 'pam60', 'pam90', 'rao', 'risler', 'structure']
# '''

# sub_matrix = get_sub_matrix(name = 'basic', gap_importance = 0)

# df = all_vs_all_pssm_cons(pssm_dict, sub_matrix= sub_matrix)                        #get the cons_scores for all pssms
# #df = all_vs_all_pssm_cons(pssm_dict,)# sub_matrix= sub_matrix)                        #get the cons_scores for all pssms
# df = add_all_vs_all(df) #compute the average for a nicer plot
# #print(df)
# # print(sub_matrix)
# # exit()
# #df2pssm_visual(df, path = './temp/temp_visual.svg', alignment_index = True, scip_gap_window = 1,drop_panalty = 0, plot_all=True, drop_class_label = True)
# df2pssm_visual(df, path = './temp/temp_visual_basic_gi_0', 
# 	alignment_index = True, 
# 	#hm_all = ['1','2'],
# 	#hm_ova = [],
# 	#hm_ava = True,
# 	#pl_all = ['1','2'],
# 	#pl_ova = ['1'],
# 	#pl_ava = True,
# 	window = 5,
# 	window_type = 'max',
# 	get_best = 20,
# 	chunksize = 100,
# 	#hm_ava = True,
# 	)
	# plot_all = True,
	# scip_gap = True, 
	# scip_gap_window = 1, 
	# drop_panalty = 0.4, 
	# sg_type = 'plot',
	# get_best = 10)

# mark_top = 10) #plot_window=[20,40])
# fig.show() 
# print(df)




# # ##############
# # # align test sequence in order to apply scoring function
# # ##############

# TEMP_ALIGNMENT_PATH = os.path.join('./temp', 'temp_alignment.fasta')
# temp_folder = './temp'
# temp_alignment = domain_extractor.discri_dict2alignment(discri_dict, returns = True)
# aligned_test_seq = domain_extractor.add_seq2alignment(single_fasta = test_seq, alignment = temp_alignment, temp_files = temp_folder, debug = False)

# scoring_dict = pssm_scoring(aligned_test_seq, pssm_dict)

# print weighed_pssm_scoring(scoring_dict, weigth_matrix = df, no_weigth_matrix = True, position_threshold = 4)
# print weighed_pssm_scoring(scoring_dict, weigth_matrix = df, no_weigth_matrix = False, position_threshold = 4)


#weighed_pssm_scoring(aligned_test_seq, pssm_dict, weigth_matrix)



# # for keys, seqs in discri_dict.iteritems():
# #     print(keys)
# #     for seq in seqs:
# #         print(seq.seq)
# #         print(len(seq.seq))

# pssm_dict = st.discri_dict2pssm_dict(discri_dict) 
# df = all_vs_all_pssm_cons(pssm_dict, sub_matrix = sub_matrix)                     #get the cons_scores for all pssms
# df = add_average(df) #compute the average for a nicer plot
# df2pssm_visual(df, path = './temp/temp_visual.svg', alignment_index = True)

# print df
# print df_sum

