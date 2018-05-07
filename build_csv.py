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

parser = argparse.ArgumentParser(
	description='''Creates the csv-file, which can be used to 
	label the classes in the sequence alignment''')

parser.add_argument(
					'-i', '--input',
					required=True,
					dest = 'input',
					help='Input file, must be a sequence alignment',
					)

parser.add_argument(
					'-a', '--alignment',
					default = 'temp_alignment.fasta',
					dest = 'ali',
					help='''Output alignment file, the input alignment
					gets converted (some extra signes must be changed),
					this file is the alignment which corresponds to the
					csv file''',
					)

parser.add_argument(
					'-r', '--regex',
					dest = 'regex',
					default = None,
					help='''
					Based on this regex the class labels can be extracted from the 
					alignement: Example: '\|([a-zA-Z0-9\-_]+)$' would extract
					from the sequence ">seq bla bla #class_A" -> "class_A".
					Should be tested with an text editor that support regex parsing,
					such as notepad or sublime
					'''
					)

parser.add_argument(
					'-o', '--output',
					dest = 'csv',
					default = 'class_labels.csv',
					help='Output csv file'
					)

parser.add_argument(
					'-d', '--delete',
					action='store_true',
					dest = 'delete',
					#default = True,
					help='Allows to overwrite the created files'
					)


args = parser.parse_args()

#execute the parsed arguments trough the pred_rule object 

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