import os
from . import basic
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def file_input_check(input_file):
	'''
	input_file -> django model file 
	In order to avoid possible errors 
	a bottleneck is introduced here.	
	Input is limited to fasta file containing pure DNA
	code (multiple entries are possbile),
	this function checks this input and tries to handle any possible error
	produced by the input
	'''

	return_list = [None,None] #error, biopython

	#check file size
	#8652835
	size = 100000000
	if input_file.size >= size:
		input_error = '''
		The submitted file is bigger then {0} bytes, 
		this is the SeMPI limit, in order to keep the server stable
		and process equal amounts of data for all users.
		If your file contains mutiple entries you can split them up and
		process them, otherwise, are you sure this is genome file ?
		Example: Streptomyces albus J1074, complete genome: 6.8 MB.
		'''.format(str(size/100000000)+' MB')
		return_list[0] = input_error 
		return(return_list)

	#leave the handling to biopython
	try:
		record = list(SeqIO.parse(input_file.path, "fasta"))
		if list(record) == []:
			return_list[0] = 'No record found, maybe a parsing error, please check the file for errors (only fasta is suported)'
		else:
			return_list[1] = record

		return(return_list)
	except Exception as bio_error:
		return_list[0] = bio_error
		return(return_list)


class Conversion_Error(Exception):
	'''custom error messages, will improve debugging'''
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

def gbk2protein_fasta(GBK_PATH, FASTA_PATH, feature_flag = 'CDS'):

	'''function which extracts all the protein translations from a GBK, for a give CDS and writes it as multiple fasta into another folder,
	very handy if specific domains should be extracted from a gene cluster in order to build new hmmer profiles'''

	basic.create_folder(FASTA_PATH, remove = True)

	for file in os.listdir(GBK_PATH):

		input_file = os.path.join(GBK_PATH, file)

		input_handle  = open(input_file, "r")
		record = SeqIO.parse(input_handle, "genbank")

		for seq_record in record:

			s_id = seq_record.id

			FASTA_FILE_PATH = os.path.join(FASTA_PATH, s_id)

			output_handle = open(FASTA_FILE_PATH, "w")

			for feature in seq_record.features:

				location = str(feature.location.start) + ' # ' + str(feature.location.end)+ ' # ' + str(feature.location.strand)
				cds_entry = ''

				if feature.type in feature_flag:

					for key in list(feature.qualifiers.keys()):

						if key == 'translation':

							for elem in feature.qualifiers[key]:

								cds_entry += elem


					output_handle.write(">"+ str(s_id) + ' # ' + location + '\n' + cds_entry + '\n')

			output_handle.close()

def intelligent_conversion(INPUT_FILE_PATH, CONVERTED_INPUT_PATH):
	'''function which tries to handle multiple possible input files, output is always fasta'''
	input_handle  = open(INPUT_FILE_PATH, "r")
	record = SeqIO.parse(input_handle, "genbank")
	if len(list(record)) != 0:
		print('genbank file found') 
		gbk2DNA_fasta(INPUT_FILE_PATH, CONVERTED_INPUT_PATH)

	input_handle  = open(INPUT_FILE_PATH, "r")
	record = SeqIO.parse(input_handle, "fasta")
	if len(list(record)) != 0:
		print('fasta file found')
		fasta2DNA_fasta(INPUT_FILE_PATH, CONVERTED_INPUT_PATH)


def gbk2DNA_fasta(INPUT_FILE_PATH, CONVERTED_INPUT_PATH):
	'''function which extracts the dna part of a gbk file, for example if a already annotated gbk file should be reannotated'''

	basic.create_folder(CONVERTED_INPUT_PATH, remove = True)

	input_handle  = open(INPUT_FILE_PATH, "r")
	record = SeqIO.parse(input_handle, "genbank")
	records = list(record)

	INPUT_FILE_NAME = INPUT_FILE_PATH.split('/')[-1]

	#this can change to the correct file ending, even if there is a unusual ending.
	file_without_ending = INPUT_FILE_NAME.split('.')[0]
	file_with_new_ending = file_without_ending + '.fasta'

	OUTPUT_FILE_NAME = os.path.join(CONVERTED_INPUT_PATH, file_with_new_ending)

	if len(records) == 0:
		raise Conversion_Error('no record could be found in the file supplyed')
	if len(records) > 1:
		raise Conversion_Error('only one GBK record as input supported')

	else:
		rec = records[0]
		SeqIO.write(rec, OUTPUT_FILE_NAME, 'fasta')

def fasta2DNA_fasta(INPUT_FILE_PATH, CONVERTED_INPUT_PATH):
	'''function which extracts the dna part of a fasta file, for example if a already annotated gbk file should be reannotated'''
	#todo: check for DNA code !

	basic.create_folder(CONVERTED_INPUT_PATH, remove = True)

	input_handle  = open(INPUT_FILE_PATH, "r")
	record = SeqIO.parse(input_handle, "fasta")
	records = list(record)

	INPUT_FILE_NAME = INPUT_FILE_PATH.split('/')[-1]

	#this can change to the correct file ending, even if there is a unusual ending.
	file_without_ending = INPUT_FILE_NAME.split('.')[0]
	file_with_new_ending = file_without_ending + '.fasta'

	OUTPUT_FILE_NAME = os.path.join(CONVERTED_INPUT_PATH, file_with_new_ending)

	if len(records) == 0:
		raise Conversion_Error('no record could be found in the file supplyed')
	if len(records) > 1:
		raise Conversion_Error('only one fasta record as input supported')

	else:
		rec = records[0]
		SeqIO.write(rec, OUTPUT_FILE_NAME, 'fasta')


def fasta_check(INPUT_FILE_PATH):
	'''checks if a fasta file is an alignment or raw dna, if it is raw dna it checks if all entries have the same size.'''

	#bad_signes = ['|', '/', '(', ')'] #some signes are not good input, so they are raplaced.
	# '(', ')' due to mafft handling

	bad_signes = [] #here no sign conversion is needed

	records = list(SeqIO.parse(INPUT_FILE_PATH, "fasta"))
	all_records = ''.join([str(sequence.seq) for sequence in records]) #a text file with all the records

	for rec in records:

		#replace bad sign
		for bad_signe in bad_signes:
			if bad_signe in rec.name:
				rec.name = rec.name.replace(bad_signe, '#')

		#rec.description = '' <- this is the complete information
		rec.id = rec.description
		rec.description = ''

	if '-' in all_records or '.' in all_records: 
		return('alignment', records)

	else:
		all_lenghts = [len(sequence.seq) for sequence in records]
		if len(set(all_lenghts)) == 1:
			return('raw', records)
		else:
			return('raw', records) #actually I think different sizes are ok, as long as the alignment is correct.
			print('records have different sizes')
