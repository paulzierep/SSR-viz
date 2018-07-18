from .basic import file_conversion
from .basic import basic
from Bio import SeqIO
from . import domain_extractor
import os

def csv_builder(input_alignment, 	#path
				output_csv,			#file
				regex_string,		#string
				output_alignment,	#file
				replace_ali,		#true
				replace_csv,		#true
				):

	print(locals())

	'''
	Creates the csv file from the alignment
	'''

	##################################################
	# Set up paths
	##################################################

	input_alignment_folder = os.path.dirname(input_alignment)
	output_csv_path = os.path.join(input_alignment_folder, output_csv)

	if replace_ali:
		output_alignment_path = input_alignment
	else:
		output_alignment_path = os.path.join(input_alignment_folder, output_alignment)

	##################################################
	# Check fasta input
	##################################################

	fasta_type, records = file_conversion.fasta_check(
	input_alignment)  # checks if raw or alignment

	if fasta_type != 'alignment':  # input is not alignment
		print('Error: Input file must be alignment')
		exit()
	else:
		records = domain_extractor.unique_sequence_list(
			records)  	# check if records are unique in the list
						# creates the alignment file
		SeqIO.write(records, output_alignment_path, 'fasta')

	##################################################
	# create csv
	##################################################

	if basic.create_file(output_csv_path, remove=replace_csv):

		if not regex_string:
			# if no regex is defined an csv with only fasta names is
			# created, which can be used to enter discri rules by hand or
			# specific function
			domain_extractor.align2csv(output_alignment_path, output_csv_path)

		else:

			discri_dict = domain_extractor.align2discri_dict(output_alignment_path, regex_str=regex_string, exclude=[])
			domain_extractor.discri_dict2csv(discri_dict, output_csv_path)

	else:
		print('csv could not be replaced, as it already exsisted')

	return(output_alignment_path, output_csv_path)
