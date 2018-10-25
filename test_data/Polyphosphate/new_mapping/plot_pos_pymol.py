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

#needs a csv with already annotated sequence
csv_files = ['SSP_plot_basic_stats.csv',
			'SSP_plot_pam30_stats.csv',
			'SSP_plot_pam300_stats.csv',
			'SSP_plot_blosum30_stats.csv',
			'SSP_plot_blosum300_stats.csv',
			'SSP_plot_pc_stats.csv',
			]


x_vs_x = 'malonyl vs methoxymalonyl'
pos_column = '2hg4.pdb'
pdb_file = pos_column
chain_id = 'A'

resi_dict = {}

for csv_f in csv_files:
	name = csv_f.replace('_stats.csv','')
	name = name.replace('SSP_plot_','')
	df_dict = read_csv(csv_f)
	df = df_dict[x_vs_x]
	print(df[pos_column])
	df = df[pos_column]
	resi_dict[name] = df.as_matrix()

print(resi_dict)





#########################
#pymol
#########################

# #!/usr/bin/env python
 
# # Tell PyMOL to launch quiet (-q), fullscreen (-e) and without internal GUI (-i)
# import __main__
# __main__.pymol_argv = [ 'pymol', '' ]
 
# import pymol
 
# # Call the function below before using any PyMOL modules.
# pymol.finish_launching()
 
# from pymol import cmd

# cmd.load(pdb_file)
# cmd.hide('all')
# cmd.show('cartoon', 'chain {0}'.format(chain_id))


# selection_str = 'chain {0} and resi {1}'.format(chain_id, resis)
# print(selection_str)

# cmd.select('aligned_area', selection_str)
# cmd.color('red','aligned_area')
# #cmd.show_as('cartoon', 'aligned_area')
# #cmd.show_as('sticks', 'aligned_area')
# cmd.zoom('aligned_area', buffer=15,state=0,complete=0)

# cmd.show('cartoon', 'aligned_area')
# cmd.show('sticks', 'aligned_area')


