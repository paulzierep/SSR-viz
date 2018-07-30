import sys
import os
import re
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
#import copy
import itertools
import copy

import csv

SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__)) #current script path.

def clustalo_alignment(INPUT_FASTA_PATH, INPUT_FASTA_ALIGNED_PATH, debug = False):
	'''creates an alignment with clustal omega'''

	sub_report = subprocess.Popen('clustalo --force --infile={0} -o {1}'.format(INPUT_FASTA_PATH, INPUT_FASTA_ALIGNED_PATH), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	
	if debug != False: #detailed reports only in debug mode
		print((sub_report[0]))
		print((sub_report[1]))
	print((INPUT_FASTA_PATH + ' aligned'))

# def clustalo_dist_matrix(INPUT_FASTA_PATH = None, INPUT_SEQ_LIST = None, MATRIX_OUTPUT = None, return_matrix = False, debug = False):

# 	TEMP_PATH_IN = './temp/raw_seqs_temp.fasta'
# 	TEMP_PATH_OUT = './temp/matrix_temp.txt'

# 	if not INPUT_FASTA_PATH and not INPUT_SEQ_LIST:
# 		print('Either fasta file or sequence list must be supplied')
# 		exit()

# 	if not MATRIX_OUTPUT and not return_matrix:
# 		print('Either matrix output or return matrix must be defined')
# 		exit()

# 	if not INPUT_FASTA_PATH:
# 		SeqIO.write(INPUT_SEQ_LIST, TEMP_PATH_IN, 'fasta')
# 		INPUT_FASTA_PATH = TEMP_PATH_IN

# 	if not MATRIX_OUTPUT:
# 		MATRIX_OUTPUT = TEMP_PATH_OUT

# 	sub_report = subprocess.Popen('clustalo --infile={0} -v --distmat-out {1}'.format(INPUT_FASTA_PATH, MATRIX_OUTPUT), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

# 	if debug != False: #detailed reports only in debug mode
# 		print(sub_report[0])
# 		print(sub_report[1])

	# if return_matrix:
	# 	return()




def clustalo_alignment_flat(INPUT_FASTA_PATH = None, INPUT_SEQ_LIST = None, OUTPUT_ALIGNED_PATH = None, return_seqs = False, debug = False):
	'''creates an alignment with clustal omega'''

	#the temp files are always created relative from this script
	TEMP_PATH_IN = os.path.join(SCRIPT_PATH, 'temp/raw_seqs_temp.fasta')
	TEMP_PATH_OUT = os.path.join(SCRIPT_PATH,'temp/align_seqs_temp.fasta')

	#set temp paths if not supplied
	if not INPUT_FASTA_PATH and not INPUT_SEQ_LIST:
		print('Either fasta file or sequence list msut be supplied')
		exit()

	if not INPUT_FASTA_PATH:
		SeqIO.write(INPUT_SEQ_LIST, TEMP_PATH_IN, 'fasta')
		INPUT_FASTA_PATH = TEMP_PATH_IN

	if  not OUTPUT_ALIGNED_PATH:
		OUTPUT_ALIGNED_PATH = TEMP_PATH_OUT

	sub_report = subprocess.Popen('clustalo --force --infile={0} -o {1}'.format(INPUT_FASTA_PATH, OUTPUT_ALIGNED_PATH), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	
	if debug != False: #detailed reports only in debug mode
		print((sub_report[0]))
		print((sub_report[1]))

	# print(INPUT_FASTA_PATH + ' aligned')

	if return_seqs:
		alignment = list(SeqIO.parse(OUTPUT_ALIGNED_PATH, 'fasta'))
		return(alignment)
	else:
		return(OUTPUT_ALIGNED_PATH)

	#k = subprocess.Popen('clustalo --force --infile={0} --outfmt=st -o {1}'.format(multi_f + '_multi_fasta', multi_f +'_align.stockholm'), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

def add_seq2alignment(	INPUT_SINGLE_FASTA_PATH = None,\
						INPUT_FASTA_ALIGNED_PATH = None,\
						single_fasta = None,\
						alignment = None,\
						OUTPUT_PATH = None,\
						temp_files = None,\
						debug = False):
	'''
	Input: 	INPUT_SINGLE_FASTA_PATH (sequence which will be added to the alignment)
			INPUT_FASTA_ALIGNED_PATH (alignment the sequence will be added)
			OUTPUT_PATH (the new alignment, with the seq included)
	
			The Biopthon_Objects of seq and alignment can also be passed, 
			they will be stored temporarely

			a temp file folder needed to be supplied, best be the pred rule folder

	Return: seq_a (the added aligned sequence)

	add a sequence to an alignment, the alignment stays the same, very good for example for a tree where the positions should not change;
	this does not work with clustalo, as clustalo can do this only via profiles, were always all inputs can be changed, therefore mafft is
	used in this case, a very nice alignemnt tool as well, the --keeplenght is used, becouse totally new insertions can not give more informations 
	due to the fact, that all the other seqs have gaps at this position

	--anysymbol is included in order to allow for AS like Selenocysteine (U)
	'''

	##########################
	# check if all data is supplied
	##########################

	if not temp_files:
		print('Temp folder needs to be supplied')
		exit()

	if INPUT_SINGLE_FASTA_PATH == None and single_fasta == None:
		print('Single fasta path or bio single fasta must be supplied.')
		exit()

	if INPUT_FASTA_ALIGNED_PATH == None and alignment == None:
		print('Alignment fasta path or bio alignment fasta must be supplied.')
		exit()

	##########################
	# create files if needed, by default the paths are overwritten if the Bio file is supplied.
	##########################

	if single_fasta:
		INPUT_SINGLE_FASTA_PATH = os.path.join(temp_files, 'input_single_fasta_temp.fasta')
		SeqIO.write(single_fasta, INPUT_SINGLE_FASTA_PATH, 'fasta')

	if alignment:
		INPUT_FASTA_ALIGNED_PATH = os.path.join(temp_files, 'input_fasta_aligned_temp.fasta')
		SeqIO.write(alignment, INPUT_FASTA_ALIGNED_PATH, 'fasta')

	if not OUTPUT_PATH:
		OUTPUT_PATH = os.path.join(temp_files, 'alignment_with_added_seq_temp.fasta')

	##########################
	# add sequence to alignment
	##########################

	sub_report = subprocess.Popen('mafft --anysymbol --addfull {0} --keeplength {1} > {2}'.format(INPUT_SINGLE_FASTA_PATH, INPUT_FASTA_ALIGNED_PATH, OUTPUT_PATH), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

	if debug != False:  #detailed reports only in debug mode
		print((sub_report[0]))
		print((sub_report[1]))

	#print(INPUT_SINGLE_FASTA_PATH + ' aligned to' + INPUT_FASTA_ALIGNED_PATH)

	single_seq = SeqIO.parse(INPUT_SINGLE_FASTA_PATH, "fasta")
	alignment = SeqIO.parse(OUTPUT_PATH, "fasta")

	##########################
	# get newly aligned sequence
	##########################

	for seq_S in single_seq:
		for seq_a in alignment:
			if seq_S.name == seq_a.name:
				return(seq_a)





###################################
# debug 
###################################

# test_seq_aligned = add_seq2alignment(INPUT_SINGLE_FASTA_PATH = './temp/test_seq.fasta',\
# 													INPUT_FASTA_ALIGNED_PATH = './temp/alignment.fasta',\
# 													OUTPUT_PATH = './temp/added_seq2alignment.fasta', \
# 													temp_files = './temp') 			#align input seq with multiple seq alignment,
# 																					#is important, so that positions are equal
# print test_seq_aligned


def alignment2raw_fasta(INPUT_FASTA_ALIGNED_PATH, INPUT_FASTA_PATH):
	'''takes an aligned seq and writes it as a raw sequence '''
	records = list(SeqIO.parse(INPUT_FASTA_ALIGNED_PATH, "fasta"))
	for sequence in records:
		sequence.seq = Seq(str(sequence.seq).replace('-','').replace('.',''))

	#print records
	SeqIO.write(records, INPUT_FASTA_PATH, 'fasta')

	print((INPUT_FASTA_PATH + ' unaligned'))


def hmm_build(INPUT_FASTA_PATH, HMMER_PATH, sub_report = False):

	'''builds a hmmer profile based on aligned sequences or single sequence'''

	sub_report = subprocess.Popen('hmmbuild --amino {0} {1}'.format(HMMER_PATH, INPUT_FASTA_PATH), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

	if sub_report == True:
		print((sub_report[0]))
		print((sub_report[1]))
		print((HMMER_PATH + ' profile build'))
	else:
		pass

def align2discri_dict(INPUT_FASTA_ALIGNED_PATH, regex_str = '_([0-9]+)$', exclude = []):
	'''
	input: 	INPUT_FASTA_ALIGNED_PATH (an aligned fasta file)
			regex_str (the regex used as base to seperate the sequence in classes)
			exclude (names to exclude from the classes)

	return: discri_dict (a dict which has the class names as keys and a list of the representing sequences as values)

	the regex defines what part in the name to use as dicrimination, numbers or names seperated with a marker at the end are easy to 
	implement...example '_([0-9]+)$' --> blabla_12 --> disriminator: 12; note that all kinds of patterns can be used, but should only occure ones in the
	name, only the part in (...) will be taken as discriminatior, only one (...) subgroup is supported'''

	align = AlignIO.read(INPUT_FASTA_ALIGNED_PATH, "fasta")
	discri_dict = {}

	for entry in align:

		#descrption is best to use, as this is the entire header in a fasta file.
		regex_search = re.search(regex_str, entry.description)

		if regex_search == None:
			print((entry.name + ' - discriminator could not be detected via regex'))
			if 'unknown' not in list(discri_dict.keys()): #create new dict entry if needed, fill list with sequences.
				discri_dict['unknown'] = []

			discri_dict['unknown'].append(entry)

		else:
			#print regex_search.group(0)
			found = regex_search.group(1)

			if found not in list(discri_dict.keys()): #create new dict entry if needed, fill list with sequences.
				discri_dict[found] = []

			discri_dict[found].append(entry)

	#print discri_dict

	for elem in exclude:
		try:
			del discri_dict[elem] #just too different.
		except:
			print((elem + ' can not be excluded'))

	#print discri_dict

	return(discri_dict)

def discri_dict2csv(discri_dict, DISCRI_CSV_PATH):
	'''dumps the discri_dict as a csv file, which can be used for further annotation'''

	with open(DISCRI_CSV_PATH, 'w') as csvfile:
		fieldnames = ['Name', 'Class']
		writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

		writer.writeheader()

		for discri in discri_dict:
			for elem in discri_dict[discri]:
				writer.writerow({'Name': elem.description, 'Class': discri})

def align2csv(INPUT_FASTA_PATH, DISCRI_CSV_PATH):
	'''if no regex is defined, the sequences can be directly transformed to csv an annotated'''

	records = list(SeqIO.parse(INPUT_FASTA_PATH, "fasta"))

	with open(DISCRI_CSV_PATH, 'w') as csvfile:
		fieldnames = ['Name', 'Class']
		writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

		writer.writeheader()

		for sequence in records:
			# print(sequence)
			# exit()
			writer.writerow({'Name': sequence.description, 'Class': ''})

# def align2discri_dict_update(INPUT_FASTA_PATH, discri_dict, regex_str = None):
# 	'''updates the discri dict with new sequences from the alignment'''

# 	records = list(SeqIO.parse(INPUT_FASTA_PATH, "fasta"))

# 	all_records = list(itertools.chain.from_iterable(list(discri_dict.values()))) #get all seqs in the discri dict
# 	all_names = [seq.name for seq in all_records] #get all the names in the discri dict

# 	for sequence in records:
# 		if sequence.name not in all_names:

# 			if regex_str == None:

# 				if 'unknown' not in list(discri_dict.keys()): #create new dict entry if needed, fill list with sequences.
# 					discri_dict['unknown'] = []

# 				discri_dict['unknown'].append(sequence)

# 			else:

# 				regex_search = re.search(regex_str, sequence.name)
# 				#print regex_search.group(1)
# 				if regex_search == None:
# 					print((sequence.name + ' - discriminator could not be detected via regex'))
# 					if 'unknown' not in list(discri_dict.keys()): #create new dict entry if needed, fill list with sequences.
# 						discri_dict['unknown'] = []

# 					discri_dict['unknown'].append(sequence)

# 				else:
# 					#print regex_search.group(0)
# 					found = regex_search.group(1)

# 					if found not in list(discri_dict.keys()): #create new dict entry if needed, fill list with sequences.
# 						discri_dict[found] = []

# 					discri_dict[found].append(sequence)				



# 	return(discri_dict)

def csv2discri_dict(DISCRI_CSV_PATH, INPUT_FASTA_ALIGNED_PATH, \
					delimiter = ',', name_field = 'Name', \
					class_field = 'Class',
					exclude=[''],
					reduce_DB=False, #must be a integer, all classes with less items are removed.
					verbose = False,
					):

	'''
	converts an exsisting csv file to the discri_dict,
	alllows to use csv_kwargs in order to pass specific parsing parameteres to the csv object
	'''

	if not os.path.exists(DISCRI_CSV_PATH):								#no csv no discri dict
		print((DISCRI_CSV_PATH ,' not found, create with get_csv'))
		return(None)

	else:

		align = list(AlignIO.read(INPUT_FASTA_ALIGNED_PATH, "fasta")) 	#read the fasta which contains the seqs
		
		discri_dict = {}

		with open(DISCRI_CSV_PATH, 'r') as csvfile:

			reader = list(csv.DictReader(csvfile, delimiter = delimiter))


			#############################
			# input check
			#############################

			if not reader:
				print('Error: Csv could not be parsed !')
				exit()

			test_reader = reader[0]

			if name_field not in list(test_reader.keys()):
				print(('Error: No matching name filed in file, maybe it is one of those: {0}'.format(list(reader.keys()))))
				exit()
			if class_field not in list(test_reader.keys()):
				print(('Error: No matching class filed in file, maybe it is one of those: {0}'.format(list(reader.keys()))))
				exit()

			#############################
			# parsing
			#############################

			if verbose:
				print('\n ***Verbose **** \n Sequence parsing: \n')

			for row in reader:

				if row[class_field] in exclude:
					if verbose:
						print((str(row[name_field]) + ' discarded, due to exclusion creteria'))
					continue

				if row[class_field] not in list(discri_dict.keys()):		#add new key to discri dict for new class
					discri_dict[row[class_field]] = []

				entry_matches = False 									#so far this entry did not match the alignment
				for entry in align:										#compare the name to the alignment
					if entry.description == row[name_field] and entry_matches == False:
						discri_dict[row[class_field]].append(entry)			#add Biopython obj if it is found in the csv
						entry_matches = True							#only one Biopython obj can be assigned to the csv entry

					elif entry.description == row[name_field] and entry_matches == True:
						if verbose:
							print(('Multiple matches between csv file and alignment' + entry.description + 'discarded'))

				if entry_matches == False:
					if verbose:
						print(('No aligned sequence found for ' + row[name_field]))
					if not discri_dict[row[class_field]]: #if the seqeuence is the only one representing this discriminator, the discri is removed from the disrci dict
						discri_dict.pop(row[class_field], None)

		if reduce_DB != False:
			pop_keys = []
			for key, val in list(discri_dict.items()):
				if len(val) < reduce_DB:
					pop_keys.append(key)
			for key in pop_keys:
				discri_dict.pop(key, None)

		return(discri_dict)

def unique_sequence_list(bio_record):
	'''checks a biopython seq list for dublicates'''
	unique_seq_list = []
	unique_rec_list = []
	for rec in bio_record:
		if rec.seq in unique_seq_list:
			print((rec.id + rec.name + rec.description + ' sequence is a dublicate and thefore discarded'))
		else:
			unique_seq_list.append(rec.seq)
			unique_rec_list.append(rec)

	#print set([len(rec) for rec in unique_rec_list])
	# 	print('sequences lenght differes')

	return(unique_rec_list)

def add_seq2raw(INPUT_FASTA_PATH, UPDATE_FASTA_PATH):
	'''updates the multi fasta file, returns True if updated, False if not (due to dublicates...)'''

	old_records = list(SeqIO.parse(INPUT_FASTA_PATH, "fasta"))
	new_records = list(SeqIO.parse(UPDATE_FASTA_PATH, "fasta"))

	#check length
	#new_records_lenght = [len(rec.seq) for rec in new_records]


	all_records = old_records + new_records #combine records
	all_records = unique_sequence_list(all_records) #unique records only

	if len(all_records) == len(old_records):
		return(False)

	else:
		SeqIO.write(all_records, INPUT_FASTA_PATH, 'fasta')
		return(True)


def discri_dict2alignment(discri_dict, OUTPUT_PATH = None, returns = False):
	'''reconverts the discri_dict to an alignment'''
	alignment = []
	for discri in discri_dict:
		for seq in discri_dict[discri]:
			alignment.append(seq)

	align1 = MultipleSeqAlignment(alignment)

	if returns == True: #the alignment can be returned and/or stored
		return(align1)

	if OUTPUT_PATH:
		SeqIO.write(align1, OUTPUT_PATH, "fasta")


def get_indices_from_sequence(seq, as_dict = {}):

	'''function which tries to find a pattern provided as dict of the form {643:'Q', 671:'L', 742:'Y', 744:'S', 795:'L'}
	in a sequence, the pattern does not have to start with zero !'''

	#############################################
	# basic input check
	#############################################

	if as_dict == {}:
		print('as_idx must be supplied !!')
		exit()

	#############################################
	# normalization, the smalest index in the pattern becomes 0, the rest is adjusted to this values
	#############################################

	min_idx = min(as_dict.keys())

	zero_idx2as_dict = {}
	#zero_as2idx_dict = {}
	for index, AS in list(as_dict.items()):
		zero_idx2as_dict[index - min_idx] = AS 

	#############################################
	# search for the sequence, if the first value is found, the other values are 
	# evaluated by creating a comparison dict (comp_dict) which is tested if equal to the providede dict.
	# if the dict matches a second dict (temp_dict) is returned, which provides the actual indices of the pattern 
	# in the sequence.
	#############################################


	index = 0
	for item in seq:
		if item == zero_idx2as_dict[0]:
			#print index

			temp_dict = {}
			comp_dict = {}

			for key in zero_idx2as_dict:
				if (index + key) <= len(seq): #take care not to excede the list
					temp_dict[index + key] = seq[index + key] #return dict
					comp_dict[key] = seq[index + key] #comaprison dict

			if comp_dict == zero_idx2as_dict:
				print('Matching sequence pattern found:')
				print(temp_dict)
				return(temp_dict)

		index += 1

	print('Matching sequence pattern could not be found')
