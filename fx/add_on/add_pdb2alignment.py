from Bio import SeqIO
import os
import shutil
import subprocess
import pandas as pd
import csv
import argparse

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import is_aa
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



#import argparse
# from gooey import Gooey
# from gooey import GooeyParser

#########################
# Parsing
#########################

# @Gooey(program_name='Add pdb')
# def main():

# 	parser = GooeyParser(

# 		description='''This programm allows to add
# pdb sequence positions to a SSR-vis csv file''')

# 	################################
# 	# CSV builder
# 	################################



# 		# metavar = 'test',
# 		# description = 'test',
# 		# help='Creates the csv-file, which can be used to label the classes in the sequence alignment')


# 	# MAFFT_PATH = '/usr/bin/mafft'

# 	# pdb_file = '1b3n.pdb'

# 	# ALIGNMENT_PATH = 'AT_alignment.fasta'

# 	# chain = 'A'

# 	# temp_folder = './tmp'

# 	# csv_input = 'SSP_plot_basic_10_stats.csv'

# 	# csv_output = 'SSP_plot_basic_10_stats_added.csv'

# 	# csv_input = 'SSP_plot_basic_10_stats_added.csv'

# 	parser.add_argument(
# 						'-a', '--ali',
# 						metavar = 'Input sequence alignment file',
# 						required=True,
# 						dest = 'ali',
# 						help='Input file must be a sequence alignment in FASTA format',
# 						widget='FileChooser'
# 						)

# 	parser.add_argument(
# 						'-p', '--pdb',
# 						metavar = 'Input pdb file',
# 						required=True,
# 						dest = 'pdb',
# 						help='Input file must be a sequence alignment in FASTA format',
# 						widget='FileChooser'
# 						)

# 	parser.add_argument(
# 						'-c', '--chain',
# 						metavar = 'Chain in the pdb file',
# 						default = 'A',
# 						required=True,
# 						dest = 'chain',
# 						help='Must exsist in the pdb',
# 						)

# 	parser.add_argument(
# 						'-t', '--temp',
# 						metavar = 'Temporary folder',
# 						required=True,
# 						# default = 'temp',
# 						dest = 'temp',
# 						help='Stores added alignment of the pdb and the mafft map file',
# 						widget="DirChooser",
# 						)

# 	parser.add_argument(
# 						'-m', '--mafft',
# 						metavar = 'Mafft executable',
# 						required=True,
# 						dest = 'mafft',
# 						default = 'mafft',
# 						help='mafft allows to add a sequence to an alignment',
			
# 						)

# 	parser.add_argument(
# 					'-i_csv', '--input_csv',
# 					metavar = 'Input csv',
# 					required=True,
# 					dest = 'i_csv',
# 					help='mafft allows to add a sequence to an alignment',
# 					widget='FileChooser'
# 					)

# 	parser.add_argument(
# 					'-o_csv', '--output_csv',
# 					metavar = 'Output csv',
# 					required=True,
# 					dest = 'o_csv',
# 					help='mafft allows to add a sequence to an alignment',
# 					)



# 	args = parser.parse_args()


# 	seq = pdb2seq(args.pdb, args.chain, args.temp)
# 	map_file = mafft_add_seq(args.mafft, args.ali, seq, args.temp)
# 	map_df = map2df(map_file)
# 	add_pos2csv(args.i_csv, args.o_csv, map_df, args.pdb)

# 	print('''
# 	#################################################################
# 	New csv file written to:
# 	{0}
# 	#################################################################
# 	'''.format(args.o_csv)
# 	)



#########################
# Add pdb functions
#########################

def mafft_add_seq(mafft_exe, alignment, seq2add, tmp):
	'''
	Add a seq to the alignment, returns a map file path
	'''
	if not shutil.which(mafft_exe):
		print('mafft executables could not be found')

	added_alignment = seq2add.replace('.fasta', '') + '_added.fasta'

	stdout,stderr = subprocess.Popen(
		'{0} --anysymbol --addfull {1} --mapout --keeplength {2} > {3}'.format(
			mafft_exe,
			seq2add,
			alignment,
			added_alignment,
			),
			shell=True,
			stdout=subprocess.PIPE,
			stderr=subprocess.PIPE
			).communicate()

	#print(stdout.decode('utf-8'))
	#print(stderr.decode('utf-8'))

	map_file = seq2add + '.map'

	return(map_file)

def map2df(map_file, seq_offset):

	with open(map_file, 'r') as os_file:
		lines = list(os_file.readlines())

	lines = lines[2:]
	lines = [line.replace('\n','') for line in lines]
	lines = [line.replace(' ','') for line in lines]
	lines = [line.split(',') for line in lines]

	labels = ['AS','pdb', 'alignment']

	df_map = pd.DataFrame.from_records(lines, columns=labels)
	#df_map['alignment'] = pd.to_numeric(df_map['alignment'], errors='coerce')
	#df_map['alignment'].fillna(-1, inplace = True)
	df_map.set_index('alignment' ,inplace = True)

	#here the pdb indices are corrected according to the offset
	df_map['pdb'] = pd.to_numeric(df_map['pdb']) + seq_offset
	return(df_map)

def ali_pos2pdb(df, positions):
	'''
	Returns the slice of the df with positions as index
	'''

	#positions = [int(pos) for pos in list(positions)]
	#print(positions)

	df = df.loc[(df.index != '-')].copy()
	df = df.reindex(positions)
	df.fillna('-', inplace = True)
	# df = df.loc[positions,:].copy()

	return(df)

def pdb2seq(pdbFile, chain, tmp):
	'''
	Returns a sequence of the chain and pdb
	'''

	# print(chain)
	# print(pdb)
	# exit()

	handle = open(pdbFile, "rU")
	
	# print(handle)
	# print(list(SeqIO.parse(handle, "pdb-seqres")))
	# exit()

	#the atom parser adds an "X" for missing indices, but 
	#needs to be adjusted to the offset.

	for record in SeqIO.parse(handle, "pdb-atom"):

		chain_pdb = record.annotations['chain'] #get the chain

		if chain == chain_pdb:

			#write the sequence in the pdb-atom
			pdb_file_name = os.path.split(pdbFile)[-1].replace('.pdb', '')
			f_name = pdb_file_name + '_' + chain
			out_file = os.path.join(tmp, f_name + '.fasta')
			with open(out_file, "w") as output_handle:
				SeqIO.write(record, output_handle, "fasta")

			#get the offset of the atom sequence

			p = PDBParser(PERMISSIVE=1)
			structure = p.get_structure(pdbFile, pdbFile)

			if (len(list(structure))) > 1:
				print('''Multiple structure parsing is not supported, plase create a pdb with only
					one structure''')
				exit()

			for model in structure:
				for pdb_chain in model:
					chainID = pdb_chain.get_id()
					if chainID == chain:
						seq_offset = ((pdb_chain.get_list()[0]).get_id()[1]) - 1

			return(out_file, seq_offset)

	print('Chain not found in the {0} pdb file'.format(pdb))

# def pdb2seq(pdbFile, chain, tmp):


# 	## First, open and parse the protein file
# 	p = PDBParser(PERMISSIVE=1)
# 	structure = p.get_structure(pdbFile, pdbFile)

# 	## Now go through the hierarchy of the PDB file
# 	##
# 	## 1- Structure
# 	##      2- Model
# 	##          3- Chains
# 	##              4- Residues
# 	##
	
# 	if (len(list(structure))) > 1:
# 		print('''Multiple structure parsing is not supported, plase create a pdb with only
# 			one structure''')
# 		exit()


# 	for model in structure:

# 		for pdb_chain in model:

# 			'''
# 			Seqres does start at the beginning of the protein, but often only a certain 
# 			part is shown in the pdb, therefore a certain number of "X" is added to the chain,
# 			so that the correct number corresponds to the sequence
# 			'''

# 			#get real starting number in the pdb
# 			seq_offset = ((pdb_chain.get_list()[0]).get_id()[1]) - 1

# 			#get gap to add to alignment
# 			#head_attachment = ''.join(['X' for x in range(number_in_pdb)])

# 			seq = list()
# 			chainID = pdb_chain.get_id()

# 			if chain == chainID:

# 				for residue in pdb_chain:
# 					## The test below checks if the amino acid
# 					## is one of the 20 standard amino acids
# 					## Some proteins have "UNK" or "XXX", or other symbols
# 					## for missing or unknown residues
# 					if is_aa(residue.get_resname(), standard=True):
# 						seq.append(three_to_one(residue.get_resname()))
# 					else:
# 						seq.append("X")
				 
# 				## This line is used to display the sequence from each chain
				 
# 				#print(">Chain_" + chainID + "\n" + str("".join(seq)))
				 
# 				## The next two lines create a sequence object/record
# 				## It might be useful if you want to do something with the sequence later
				
# 				final_seq = (str(''.join(seq)))
# 				#final_seq = #head_attachment + (str(''.join(seq)))

# 				myProt = Seq(final_seq, IUPAC.protein)
# 				record = SeqRecord(myProt, id=chainID, name="", description="")
				
# 				pdb_file_name = os.path.split(pdbFile)[-1].replace('.pdb', '')
# 				f_name = pdb_file_name + '_' + chain

# 				out_file = os.path.join(tmp, f_name + '.fasta')
# 				with open(out_file, "w") as output_handle:
# 					SeqIO.write(record, output_handle, "fasta")

# 				return(out_file, seq_offset)

# 		print('Chain not found in the {0} pdb file'.format(pdbFile))


def add_pos2csv(csv_input, csv_output, map_df, pdb_path):
	pdb_folder_ , pdb_file = os.path.split(pdb_path)

	with open(csv_input) as csvfile:
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

	pandas_dict = {}
	for key, row in row_dict.items():
		df = pd.DataFrame.from_records(row, columns = header_dict[key])
		df.set_index('Position', inplace = True)
		df.sort_index(inplace=True, ascending = False)
		pandas_dict[key] = df

	csv_text = ''
	for clas, df in pandas_dict.items():
		positions = df.index

		selection = ali_pos2pdb(map_df, positions)

		# print(selection['pdb'].index)
		# print(df.index)

		df[pdb_file.replace('.pdb', '')] = selection['pdb']
		# print(selection)
		# print(df)
		# exit()
		csv_text += clas + '\n'
		csv_text += df.to_csv()

		with open(csv_output, 'w') as csvfile:
			csvfile.write(csv_text)


# if __name__ == "__main__":
# 	main()

########################
# debug
########################

# MAFFT_PATH = '/usr/bin/mafft'

# pdb_file = '1b3n.pdb'

# ALIGNMENT_PATH = 'AT_alignment.fasta'

# chain = 'A'

# temp_folder = './tmp'

# csv_input = 'SSP_plot_basic_10_stats.csv'

# csv_output = 'SSP_plot_basic_10_stats_added.csv'

# csv_input = 'SSP_plot_basic_10_stats_added.csv'

# seq = pdb2seq(pdb_file, chain, temp_folder)
# map_file = mafft_add_seq(MAFFT_PATH, ALIGNMENT_PATH, seq, temp_folder)
# map_df = map2df(map_file)
# add_pos2csv(csv_input, csv_output, map_df, pdb_file)




