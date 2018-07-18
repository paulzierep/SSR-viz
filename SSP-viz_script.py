#!/usr/bin/env python3


import os
import argparse
#from pssm_lib.main_pred import pred_rule

#for deployment the default warnings which are
#due to biopython and notebook are switched off

import sys
import warnings

if not sys.warnoptions:
	warnings.simplefilter("ignore")

#import the prediction_rule object and give some 
#arguments via commandline

from prediction_rules.main_pred import pred_rule

from add_on import add_pdb2alignment

from gooey import Gooey
from gooey import GooeyParser

from subprocess import Popen, PIPE

#import tempfile
import json

@Gooey(program_name='SSR-viz')
def main():

	parser = GooeyParser(

		description='''The CSV_builder protocol can be used to create the 
csv annotation file, the SSR_plot protocol starts 
the plotting window''')

	subs = parser.add_subparsers(help='commands', dest='command')

	################################
	# CSV builder
	################################

	csv_parser = subs.add_parser('CSV_Builder')
		# metavar = 'test',
		# description = 'test',
		# help='Creates the csv-file, which can be used to label the classes in the sequence alignment')

	csv_parser_r = csv_parser.add_argument_group('Required argument')

	csv_parser_r.add_argument(
						'-i', '--input',
						metavar = 'Input sequence alignment file',
						required=True,
						dest = 'input',
						help='Input file must be a sequence alignment in FASTA format',
						widget='FileChooser'
						)

	csv_parse_o = csv_parser.add_argument_group('Optional arguments')

	csv_parse_o.add_argument(
						'-ca', '--convert_alignment',
						metavar = 'Inplace FASTA conversion',
						action='store_false',
						default = True,
						dest = 'ca',
						help='''Removes duplicates from the input FASTA file, does not create a temporary alignment''',
						)

	csv_parse_o.add_argument(
						'-a', '--alignment',
						metavar = 'Temporary alignment file name',
						default = 'temp_alignment.fasta',
						dest = 'ali',
						help='''Temporary alignment file without duplicates''',
						)

	csv_parse_o.add_argument(
						'-r', '--regex',
						metavar = 'Regex extraction of the class label',
						dest = 'regex',
						default = None,
						help='''Based on this regex the class labels can be extracted from the 
alignment sequence names: 
Example: >sequence1 class: I
Regex: class: ([a-z]*)$

Example2: >sequence1 class: 5
Regex2: class: ([0-9]*)$

Examples can be found in the Manuel, Section Regex.
Can be tested with an text editor that support regex parsing,
such as notepad or sublime or online: https://regex101.com/
						'''
						)

	csv_parse_o.add_argument(
						'-o', '--output',
						metavar = 'Output file for the class labels',
						dest = 'csv',
						default = 'class_labels.csv',
						help='The name of the CSV file, this is used for the SSR algorithm.'
						)

	csv_parse_o.add_argument(
						'-d', '--delete',
						metavar = 'Overwrite existing files',
						action='store_true',
						dest = 'delete',
						#default = True,
						help='Allows to overwrite the created files'
						)

	################################
	# SSR plot parser
	################################

	ssp_plot_parser = subs.add_parser('SSR_plot', help='Creates the differnce pssm plot, based on the class_label.csv file and the converted alignmet file')

	ssp_plot_parser_r = ssp_plot_parser.add_argument_group('Required arguments')

	ssp_plot_parser_r.add_argument(
						'-i', '--input-csv',
						metavar = 'Input CSV file',
						required=True,
						dest = 'csv',
						help='Must be a csv file that with corresponding names to the alignment',
						widget='FileChooser',
						)

	ssp_plot_parser_r.add_argument(
						'-cl', '--class_label',
						metavar = 'Class label',
						#required=True,
						default = 'Class',
						dest = 'cl',
						help='An alternative column with the class label in the CSV file can be specified',
						)

	ssp_plot_parser_r.add_argument(
						'-a', '--alignment',
						metavar = 'Alignment',
						required=True,
						dest = 'ali',
						help='Input file, must be a sequence alignment with corresponding names to the csv file',
						widget='FileChooser'
						)


	################################
	# SSR add pdb to MSA
	################################

	add_pdb_parser = subs.add_parser('Add_pdb', help='Add the indices of an pdb file to the created csv files')

	add_pdb_parser_r = add_pdb_parser.add_argument_group('Required arguments')

	add_pdb_parser_r.add_argument(
						'-a', '--ali',
						metavar = 'Input sequence alignment file',
						required=True,
						dest = 'ali',
						help='Input file must be a sequence alignment in FASTA format',
						widget='FileChooser'
						)

	add_pdb_parser_r.add_argument(
						'-p', '--pdb',
						metavar = 'Input pdb file',
						required=True,
						dest = 'pdb',
						help='Input file must be a sequence alignment in FASTA format',
						widget='FileChooser'
						)

	add_pdb_parser_r.add_argument(
						'-c', '--chain',
						metavar = 'Chain in the pdb file',
						default = 'A',
						required=True,
						dest = 'chain',
						help='Must exsist in the pdb',
						)

	add_pdb_parser_r.add_argument(
						'-t', '--temp',
						metavar = 'Temporary folder',
						required=True,
						# default = 'temp',
						dest = 'temp',
						help='Stores added alignment of the pdb and the mafft map file',
						widget="DirChooser",
						)

	add_pdb_parser_r.add_argument(
						'-m', '--mafft',
						metavar = 'Mafft executable',
						required=True,
						dest = 'mafft',
						default = 'mafft',
						help='mafft allows to add a sequence to an alignment',
			
						)

	add_pdb_parser_r.add_argument(
					'-i_csv', '--input_csv',
					metavar = 'Input csv',
					required=True,
					dest = 'i_csv',
					help='mafft allows to add a sequence to an alignment',
					widget='FileChooser'
					)

	add_pdb_parser_r.add_argument(
					'-o_csv', '--output_csv',
					metavar = 'Output csv',
					required=True,
					dest = 'o_csv',
					help='mafft allows to add a sequence to an alignment',
					)

	args = parser.parse_args()

	##################################
	# Call the csv building routine
	##################################

	if args.command == 'CSV_Builder':

		ALIGNMENT_PATH = os.path.join(os.path.dirname(args.input), args.ali)
		CSV_PATH = os.path.join(os.path.dirname(args.input), args.csv)

		file_path_dict = 	{	
							'multi_fasta':None, 
							'pred_eval':None,
							'ref_seq':None,
							'hmmer': None,
							'pred_rules':None,
							'alignment': ALIGNMENT_PATH, 
							'discri_csv':CSV_PATH, 
							}

		#print(args.ca)
		if not args.ca:

			print('''
			#################################################################
			Alignment file is being converted, output file:
			{0}
			#################################################################
			'''.format(args.input)
			)


			pr = pred_rule(file_path_dict, args.input, remove_output = args.delete, replace_ali = True) #init the object

		else:
			print('''
			#################################################################
			Alignment file is being converted, output file:
			{0}
			#################################################################
			'''.format(ALIGNMENT_PATH)
			)

			pr = pred_rule(file_path_dict, args.input, remove_output = args.delete) #init the object

		print('''
			#################################################################
			Class label (csv) file is being created, output file:
			{0}
			#################################################################
			'''.format(CSV_PATH)
			)

		pr.get_csv(descri_regex = args.regex, remove = args.delete) #create a csv file

	##################################
	# Call the ssp plot routine
	##################################

	if args.command == 'SSR_plot':

		OUTPUT_PATH = os.path.join(os.path.dirname(args.csv))#, args.pred_rules)

		file_path_dict = 	{	
							'multi_fasta':None, 
							'pred_eval':None,
							'ref_seq':None,
							'hmmer': None,
							'pred_rules':None,
							'alignment': args.ali, 
							'discri_csv':args.csv, 
							}


		print('''
			#################################################################
			Class label (csv) file
			-> {0}
			is being synchronised with the alignment file:
			-> {1}
			based on the csv column:
			-> {2}
			#################################################################
			'''.format(args.csv, args.ali, args.cl)
			)


		pr = pred_rule(file_path_dict, args.ali, remove_output = False) #init the object
		pr.update_descri_dict(delimiter=',', name_field='Name', class_field=args.cl, exclude=[''])

		########################################
		# call the next gooey as subprocess from here 
		# should work on any system
		########################################

		default_args = {}
		default_args['alignment'] = args.ali
		default_args['class'] = args.cl
		default_args['csv'] = args.csv
		default_args['classes'] = list(pr.discri_dict.keys())


				#write a temp file with the default args as json
		CURRENT_PATH = os.path.dirname(sys.argv[0])
		TEMP_PATH = os.path.join(CURRENT_PATH, 'ssp_viz_temp_params.txt')
		with open(TEMP_PATH, 'w') as tmp:
			json.dump(default_args, tmp)

		#execute the next gooey window in seperate process
		#if the files are scripts this should be the python script
		#if executeables this should be the exe !
		# process_call = os.path.join(CURRENT_PATH, 'SSR-viz-draw')
		# process = Popen(process_call, stdout=PIPE, stderr=PIPE)
		PYTHON_PATH = sys.executable
		# print(PYTHON_PATH)
		# print('hi')
		process = Popen([PYTHON_PATH, 'SSR-viz-draw.py'], stdout=PIPE, stderr=PIPE)
		output, error = process.communicate()
		# print(output)
		# print(error)

	##################################
	# Call add pdb routine
	##################################

	if args.command == 'Add_pdb':

		seq = add_pdb2alignment.pdb2seq(args.pdb, args.chain, args.temp)
		map_file = add_pdb2alignment.mafft_add_seq(args.mafft, args.ali, seq, args.temp)
		map_df = add_pdb2alignment.map2df(map_file)
		add_pdb2alignment.add_pos2csv(args.i_csv, args.o_csv, map_df, args.pdb)

		print('''
		#################################################################
		New csv file written to:
		{0}
		#################################################################
		'''.format(args.o_csv)
		)


if __name__ == "__main__":
	main()
