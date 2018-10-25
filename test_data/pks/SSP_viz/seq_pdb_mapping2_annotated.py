import os
import subprocess
import pandas as pd

from Bio import SeqIO

from biopandas.pdb import PandasPdb

def map_pos2csv(MAIN_FOLDER,INPUT_ALI, positions, INPUT_PDB, chain_id):

	pdb_file = os.path.join(MAIN_FOLDER, INPUT_PDB)
	ali_file = os.path.join(MAIN_FOLDER, INPUT_ALI)

	output_seq_file = INPUT_PDB.replace('.pdb', '_seq.fasta')
	output_seq_file = os.path.join(MAIN_FOLDER, output_seq_file)

	added_seq_alignment_file = INPUT_ALI.replace('.fasta', '') + '_' + INPUT_PDB.replace('.pdb','') + '.fasta'
	added_seq_alignment_file = os.path.join(MAIN_FOLDER, added_seq_alignment_file)

	mapoutput = output_seq_file + '.map' #thats how mafft does it !

	pdb_name = INPUT_PDB.replace('.pdb','')

	# Initialize a new PandasPdb object
	# and fetch the PDB file from rcsb.org
	ppdb = PandasPdb()
	ppdb.read_pdb(pdb_file)
	seq_df = ppdb.amino3to1()
	seq_df['residue_number'] = ppdb.df['ATOM']['residue_number'] #map the sequence AS to the residue number


	chain_seq = seq_df.loc[(seq_df.chain_id == chain_id)]
	chain_seq_string = ''.join(list(chain_seq.residue_name))

	pdb_seq_file = '>{0}\n{1}'.format(pdb_name, chain_seq_string)


	if not os.path.exists(output_seq_file):

		with open(output_seq_file, 'w') as os_file:
			os_file.write(pdb_seq_file)

	if not os.path.exists(mapoutput) or not os.path.exists(added_seq_alignment_file):

		stdout,stderr = subprocess.Popen(
			'mafft --anysymbol --addfull {0} --mapout --keeplength {1} > {2}'.format(
				output_seq_file,
				ali_file,
				added_seq_alignment_file),
				shell=True,
				stdout=subprocess.PIPE,
				stderr=subprocess.PIPE
				).communicate()

	# print(stdout)
	# print(stderr)

	with open(mapoutput, 'r') as os_file:
		lines = list(os_file.readlines())

	lines = lines[2:]
	lines = [line.replace('\n','') for line in lines]
	lines = [line.replace(' ','') for line in lines]
	lines = [line.split(',') for line in lines]

	labels = ['AS','pdb_seq', 'new_seq']

	df_map = pd.DataFrame.from_records(lines, columns=labels)
	chain_seq = chain_seq.reset_index()

	mapped_alignment = pd.concat([df_map, chain_seq], axis = 1)
	ali_res = mapped_alignment.loc[(mapped_alignment.new_seq != '-')]
	ali_res = ali_res.residue_number


	# print('++++++++++++++++')
	# print(positions)
	# print('++++++++++++++++')

	selection_in_ali = positions
	selection_in_ali = [str(num) for num in selection_in_ali]

	# print('++++++++++++++++')
	# print(mapped_alignment)
	# print('++++++++++++++++')

	selection = mapped_alignment.loc[(mapped_alignment.new_seq.isin(selection_in_ali))]

	# print('++++++++++++++++')
	# print(selection)
	# print(selection.residue_number)
	# print('++++++++++++++++')

	return(selection)


MAIN_FOLDER = '.'

INPUT_ALI = 'AT_alignment.fasta'

chain_id = 'A'

#csv_output = '3czp 1 vs 2'

pdb_list = ['2hg4.pdb']#'6apg.pdb']

import csv

csv_files = ['SSP_plot_basic_stats.csv',
			'SSP_plot_pam30_stats.csv',
			'SSP_plot_pam300_stats.csv',
			'SSP_plot_blosum30_stats.csv',
			'SSP_plot_blosum100_stats.csv',
			'SSP_plot_pc_stats.csv',
			]

for csv_file in csv_files:

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

	csv = ''
	for clas, df in pandas_dict.items():
		for pdb in pdb_list:
			positions = df.index
			selection = map_pos2csv(MAIN_FOLDER,INPUT_ALI, positions, pdb, chain_id)
			selection.set_index('new_seq', inplace = True)
			selection.sort_index(inplace=True, ascending = False)
			df[pdb] = selection['residue_number']
		csv += clas + '\n'
		csv += df.to_csv()

		with open(csv_file, 'w') as csvfile:
			csvfile.write(csv)