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

from gooey import Gooey
from gooey import GooeyParser

from subprocess import Popen, PIPE

import tempfile
import json

@Gooey
def main():

	parser = GooeyParser(
		description='''Creates the differnce pssm plot, based on the
					class_labela.csv file and the converted alignmet file''')


	parser.add_argument(
						'-i', '--input-csv',
						required=True,
						dest = 'csv',
						help='Input file, must be a csv file that with corresponding names to the alignment',
						widget='FileChooser',
						)

	parser.add_argument(
						'-cl', '--class_label',
						#required=True,
						default = 'Class',
						dest = 'cl',
						help='An alternative column with the class label in the csv can be specified',
						)

	parser.add_argument(
						'-a', '--alignment',
						required=True,
						dest = 'ali',
						help='Input file, must be a sequence alignment with corresponding names to the csv file',
						widget='FileChooser'
						)


	args = parser.parse_args()

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
	TEMP_PATH = os.path.join('temp', 'temp_class.txt')
	with open(TEMP_PATH, 'w') as tmp:
		json.dump(default_args, tmp)

	#execute the next gooey window in seperate process
	PYTHON_PATH = sys.executable
	process = Popen([PYTHON_PATH, 'plot_dpssm_draw.py'], stdout=PIPE, stderr=PIPE)
	output, error = process.communicate()
	print(output)
	print(error)


if __name__ == "__main__":
	main()
