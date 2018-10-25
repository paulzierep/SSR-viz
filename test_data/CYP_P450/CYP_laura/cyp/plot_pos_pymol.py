import csv
import pandas as pd

def read_csv(csv_file):
	with open(csv_file) as csvfile:
		row_dict = {}
		header_dict = {}
		reader = csv.reader(csvfile)
		for row in reader:
			if 'vs' in row[0]:
				vs_key = row[0]
				row_dict[row[0]] = []
			elif 'Position' in row[0]:
				columns = row
			else:
				header_dict[vs_key] = columns
				row_dict[vs_key].append(row)

	# print(header_dict)

	pandas_dict = {}
	for key, row in row_dict.items():
		df = pd.DataFrame.from_records(row, columns = header_dict[key])
		df.set_index('Position', inplace = True)
		df.sort_index(inplace=True, ascending = False)
		pandas_dict[key] = df

	return(pandas_dict)

##############
# different matrix
#############

# #needs a csv with already annotated sequence
# csv_files = ['SSP_plot_basic_stats.csv',
# 			'SSP_plot_pam30_stats.csv',
# 			'SSP_plot_pam300_stats.csv',
# 			'SSP_plot_blosum30_stats.csv',
# 			'SSP_plot_blosum100_stats.csv',
# 			'SSP_plot_pc_stats.csv',
# 			]


# x_vs_x = 'malonyl vs methoxymalonyl'
# pos_column = '2hg4.pdb'
# pdb_file = pos_column
# chain_id = 'A'

# resi_dict = {}

# for csv_f in csv_files:
# 	name = csv_f.replace('_stats.csv','')
# 	name = name.replace('SSP_plot_','')
# 	df_dict = read_csv(csv_f)
# 	df = df_dict[x_vs_x]
# 	#print(name)
# 	df = pd.to_numeric(df[pos_column]).dropna()
# 	resi_dict[name] = df.as_matrix()

##############
# x_vs_x
#############

# #needs a csv with already annotated sequence
csv_files = [
			'cyp_basic_gi0_stats.csv',
			]


#x_vs_x = 'blue vs green'
pos_column = 'bicC_model1_robetta.pdb'
pdb_file = pos_column
chain_id = 'A'

resi_dict = {}

for csv_f in csv_files:
	name = csv_f.replace('_stats.csv','')
	name = name.replace('SSP_plot_','')
	df_dict = read_csv(csv_f)
	#df = df_dict[x_vs_x]
	for name, df in df_dict.items():
	#print(name)
		df = pd.to_numeric(df[pos_column]).dropna()
		resi_dict[name] = df.as_matrix()

##############
# to latex
#############

#needs a csv with already annotated sequence
# csv_files = [
# 			'SSP_plot_basic_10_stats.csv',
# 			]


# x_vs_x = 'malonyl vs methoxymalonyl'
# pos_column = '2hg4.pdb'
# pdb_file = pos_column
# chain_id = 'A'

# resi_dict = {}

# for csv_f in csv_files:
# 	name = csv_f.replace('_stats.csv','')
# 	name = name.replace('SSP_plot_','')
# 	df_dict = read_csv(csv_f)
# 	#df = df_dict[x_vs_x]
# 	for name, df in df_dict.items():
# 		#print(df)
# 		# exit()
# 		# df = df.iloc[[1,2,4,6]]
# 		# print(df)
# 		# exit()
# 		#df = df[['Score', pos_column, 'ethylmalonyl Cons',  ]].copy()
# 		df.sort_values(by=['Score'], inplace = True, ascending = False)
# 		df.set_index(pos_column, inplace = True)
# 		df['Score'] = df['Score'].astype('float').round(3)
# 		print(name)
# 		print(df.to_latex())
# 		#resi_dict[name] = df.as_matrix()

# exit()

#########################
#pymol
#########################

#!/usr/bin/env python
 
# Tell PyMOL to launch quiet (-q), fullscreen (-e) and without internal GUI (-i)
import __main__
__main__.pymol_argv = [ 'pymol', '' ]
 
import pymol
 
# Call the function below before using any PyMOL modules.
pymol.finish_launching()
 
from pymol import cmd

cmd.load(pdb_file)
cmd.hide('all')
cmd.show('cartoon', 'chain {0}'.format(chain_id))

for name, resis in resi_dict.items():

	resis = list(resis.astype('int').astype('str'))
	#resis = ['' if r == '' else r for r in resis]
	resis = '+'.join(resis)

	#print(resis)
	selection_str = 'chain {0} and resi {1}'.format(chain_id, resis)
	#print(selection_str)

	cmd.select(name, selection_str)
	#cmd.create(name + '_copy', selection_str)#, string selection
	#cmd.label(name + '_copy', 'resi')

# cmd.color('red','aligned_area')

#cmd.show_as('cartoon', 'aligned_area')
#cmd.show_as('sticks', 'aligned_area')
#cmd.zoom('aligned_area', buffer=15,state=0,complete=0)

#cmd.show('cartoon', 'aligned_area')
#cmd.show('sticks', 'aligned_area')


