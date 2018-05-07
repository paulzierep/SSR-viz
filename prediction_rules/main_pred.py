'''
Entry point into the SeMPI classification API
@author: Paul Zierep
'''


import os
from Bio import SeqIO

import pickle

import copy

# custom functions
from . import domain_extractor

import sys
#sys.path.append('../')
from .basic import basic
from .basic import file_conversion

from . import prediction_rules
#import eval_pred_rule
#import reference

#from multiclass_ROC import multi_rule_AUC

from collections import OrderedDict

# natural sorting makes the order of names look much nice !
from natsort import natsorted

class pred_rule():
	'''class wich stores the file and file_paths of the main files needed 
	for creation of prediction rules'''

	def __init__(self, file_path_dict, init_file, remove_output=False):

		self.allowed_paths = ['multi_fasta',
							  'alignment',
							  'hmmer',
							  'discri_csv',
							  'pred_rules',
							  'pred_eval',
							  'ref_seq',
							  'descri_dict_output']

		self.prediction_rules = {}  # name:prediction rule object
		self.prediction_rules_eval = {}  # name:prediction rule eval object
		self.reference_seqs = {}  # not jet impplemented

		#############################
		# the instance is bound to certain files, if the update function is 
		#called these files are updated
		#############################

		for key, value in list(file_path_dict.items()):
			if key not in self.allowed_paths:  # check if the path is allowed
				print((key + ' not in allowed paths'))

			else:
				# create a class variable for each path
				setattr(self, key, value)

		self.init_file_path = init_file  # get the input file

		#############################
		# initiate multi alignment and multi fasta
		#############################

		if self.multi_fasta:

		###########################
		# If not specified the input type is guessed
		###########################

			try:

				fasta_type, records = file_conversion.fasta_check(
					self.init_file_path)  # checks if raw or alignment

				if fasta_type == 'alignment':  # input is alignment

					# check if file does not already exsist
					if basic.create_file(self.alignment, remove=remove_output):
						records = domain_extractor.unique_sequence_list(
							records)  # check if records are unique in the list
						# creates the alignment file
						SeqIO.write(records, self.alignment, 'fasta')
					if basic.create_file(self.multi_fasta, remove=remove_output):
						domain_extractor.alignment2raw_fasta(
							self.alignment, self.multi_fasta)  # creates the raw file

				else:  # input is raw fasta

					# check if file does not already exsist
					if basic.create_file(self.multi_fasta, remove=remove_output):
						records = domain_extractor.unique_sequence_list(
							records)  # check if records are unique in the list
						# creates the raw file
						SeqIO.write(records, self.multi_fasta, 'fasta')
					# check if file does not already exsist
					if basic.create_file(self.alignment, remove=remove_output):
						domain_extractor.clustalo_alignment(
							self.multi_fasta, self.alignment)  # creates the alignment file

			except Exception as error:
				print(error)
				print('init_file could not be parsed')


		else:

			###########################
			# If no multi-fasta only alignment is taken
			###########################

			fasta_type, records = file_conversion.fasta_check(
			self.init_file_path)  # checks if raw or alignment

			if fasta_type != 'alignment':  # input is alignment
				print('If multi fasta is not defined input file must be alignment')
			else:
				if basic.create_file(self.alignment, remove=remove_output):
					records = domain_extractor.unique_sequence_list(
						records)  # check if records are unique in the list
					# creates the alignment file
					SeqIO.write(records, self.alignment, 'fasta')

		#############################
		# initiate hmmer profile
		#############################

		if self.hmmer:
			try:

				if basic.create_file(self.hmmer, remove=remove_output):
					domain_extractor.hmm_build(self.alignment, self.hmmer)

			except Exception as error:
				print(error)
				print('hmmer could not be build')

		#############################
		# initiate prediction rules
		#############################


		if self.pred_rules:
			basic.create_folder(
				self.pred_rules, remove=remove_output, new_folder=False)

		if self.pred_eval:
			basic.create_folder(
				self.pred_eval, remove=remove_output, new_folder=False)

		if self.ref_seq:
			basic.create_folder(
				self.ref_seq, remove=remove_output, new_folder=False)

	def get_csv(self, descri_regex, remove=False):
		'''
		initiate csv_file, either empty to create the annotation, 
		or automatic with regex discrimination
		
		the file can overwrite the old csv file, but batter 
		call it with remove == False, so that an already 
		created csv file is not removed
		'''

		self.descri_regex = descri_regex

		bad_signes = ['\|', '\/', '\(', '\)'] 

		#some signes are not good input, so they are raplaced.
		#based on the bad_signes in the file conversion
		# '(', ')' due to mafft handling

		if self.descri_regex != None:
			for bad_signe in bad_signes:
				self.descri_regex = self.descri_regex.replace(bad_signe, '\#')

		#print(self.descri_regex)

		if basic.create_file(self.discri_csv, remove=remove):

			if self.descri_regex == None:
				# if no regex is defined an csv with only fasta names is
				# created, which can be used to enter discri rules by hand or
				# specific function
				self.discri_dict = None
				domain_extractor.align2csv(self.alignment, self.discri_csv)

			else:

				try:
					'''todo exclude option'''

					self.discri_dict = domain_extractor.align2discri_dict(
						self.alignment, regex_str=self.descri_regex, exclude=[])
					domain_extractor.discri_dict2csv(
						self.discri_dict, self.discri_csv)

				except Exception as error:
					print(error)
					print('csv could not be build')

		else:
			print('csv could not be replaced, as it already exsisted')

	# def include_sequences(self, update_sequences_path):

	# 	'''
	# 	new_sequences can be included, they will be added to all the files
	# 	'''

	# 	# new seqs are included into the multi faster
	# 	update_status = domain_extractor.add_seq2raw(
	# 		self.multi_fasta, update_sequences_path)

	# 	if update_status == True:  # alignments and hmmer are only build if there are indeed new sequences

	# 		domain_extractor.clustalo_alignment(
	# 			self.multi_fasta, self.alignment)
	# 		domain_extractor.hmm_build(self.alignment, self.hmmer)

	def update_descri_dict(self, delimiter=',', name_field='Name', class_field='Class', exclude=[''], reduce_DB=False):

		'''
		updates the descri_dict, based on the annotation in the discri_csv file
		'''

		self.discri_dict = domain_extractor.csv2discri_dict(self.discri_csv,
															self.alignment,
															delimiter=delimiter,
															name_field=name_field,
															class_field=class_field,
															exclude=exclude,
															reduce_DB=reduce_DB,)  # update discri dict with current annotation


		# print(self.discri_dict)
		# exit()
	# def update_csv(self, descri_regex=None):

	# 	'''
	# 	takes the current csv, appends aligned seqs which are not yet included.
	# 	(best to be called after include_sequences)
	# 	'''

	# 	self.discri_dict = domain_extractor.csv2discri_dict(
	# 		self.discri_csv, self.alignment)  # update discri dict with current annotation

	# 	# if no regex is defined, the csv is updated with only the name of the
	# 	# new sequences
	# 	self.discri_dict = domain_extractor.align2discri_dict_update(
	# 		self.alignment, self.discri_dict, regex_str=descri_regex)
	# 	domain_extractor.discri_dict2csv(self.discri_dict, self.discri_csv)

	def update_prediction_rule(self, prediction_rule, prediction_name = None, **kwargs):

		'''
		create a specific prediction rule
		'''

		# if the name is not specified it is automatically the rule type
		if not prediction_name:
			prediction_name = prediction_rule

		if self.discri_dict == None:
			print('Discri_dict not jet computed')
		else:
			rule_class = prediction_rules.prediction_dispatcher(
				prediction_rule)  # get the rule based on the name
			rule_object = rule_class(prediction_name, self.pred_rules)
			rule_object.setup_prediction_data(**kwargs)
			rule_object.create_prediction_data(self.discri_dict, **kwargs)

			if not basic.check_dict_key_true(kwargs, 'no_rule'):
				pickle.dump(rule_object, open(rule_object.rule_storage, "wb"))

	def load_prediction_rules(self):

		''''
		load all the computed prediction rules
		'''

		# rule_path = os.path.join(self.pred_rules, prediction_rule + '.rule')

		if len(os.listdir(self.pred_rules)) == 0:
			print(('no prediction rules computed in ' + self.prediction_rules))
		else:

			for file in os.listdir(self.pred_rules):
				if '.rule' in file:
					rule_path = os.path.join(self.pred_rules, file)
					rule = pickle.load(open(rule_path, "rb"))
					self.prediction_rules[file.replace('.rule', '')] = rule

	def update_prediction_rule_evaluation(self, 
		prediction_rule,
		eval_type='loo',
		rule_name = None,
		**kwargs):

		'''
		create the data and eval object which is needed to score the 
		rule object

		the loo evaluation method only works if more then one sequence
		are representative for each discriminator
		'''

		# if the name is not specified it is automatically the rule type
		if not rule_name:
			rule_name = prediction_rule

		if self.discri_dict == None:
			print('Discri_dict not jet computed')

		else:
			evaluation_object = eval_pred_rule.eval_rule(
				prediction_rule,
				rule_name,
				eval_type,
				self.pred_eval,
				self.discri_dict, 
				**kwargs)

			evaluation_object.create_eval_data()
			print((evaluation_object.storage + ' data created'))
			pickle.dump(evaluation_object, open(
				evaluation_object.storage + '.eval', "wb"))
			print((evaluation_object.storage + '.eval created'))

	def load_prediction_rule_evaluations(self):

		'''
		loads all evaluation object, can be used for score comparison
		'''

		if len(os.listdir(self.pred_eval)) == 0:
			print(('no prediction rules computed in ' + self.pred_eval))

		else:

			for file in os.listdir(self.pred_eval):
				if '.eval' in file:
					eval_path = os.path.join(self.pred_eval, file)
					evaluation = pickle.load(open(eval_path, "rb"))
					self.prediction_rules_eval[
						file.replace('.eval', '')] = evaluation

	def compare_prediction_rules(self, rules_list=[], mode='micro', auc_name = 'temp_AUC', **AUC_kwargs):

		'''
		compares all specified rules with an ROC cuve and the specified mode
		'''

		# if rules_list == []:
		# 	print('rules must be selected')

		#ordered dict keeps the order of the keys as they are inserted,
		#good for the plot titles
		fpr_tpr_dict = OrderedDict()
		titles = []

		rules_list = natsorted(rules_list) #order by alphabet and number (1,5,10) not (1,10,5)

		for pred_rule in rules_list:
			if pred_rule in self.prediction_rules_eval:

		# for pred_rule in self.prediction_rules_eval:
		#     if pred_rule in rules_list:
				if mode == 'micro':
					fpr_tpr_dict[pred_rule] = self.prediction_rules_eval[
						pred_rule].get_AUC(return_micro=True)
				if mode == 'macro':
					fpr_tpr_dict[pred_rule] = self.prediction_rules_eval[
						pred_rule].get_AUC(return_macro=True)
				#titles.append(pred_rule)

			else:
				print(('{0} not found in the loaded prediction_rules'.format(pred_rule)))


		title = '{0} averaged ROC curves of each algorithm'.format(
			mode.capitalize())

		save_path = os.path.join(self.pred_eval, '{0}.pdf'.format(auc_name))
		multi_rule_AUC(fpr_tpr_dict, title, save_path, **AUC_kwargs)










	#############################
	# here reference sequences can be added, which can be used for example for a position
	# specific pssm, not jet tested and implemented correctly
	#############################

	def update_reference_seqs(self, name, reference_seq, file_type):
		'''
		name is the working name used in the dict., reference_seq where the initial file can be found,
		file_type is either fasta or pdb
		'''

		#ref_path = os.path.join(self.ref_seq , name + '.pickle')
		ref_object = reference.reference_sequence(
			name, self.ref_seq, reference_seq, file_type)
		ref_object.get_aligned_chain(self.alignment)

		pickle.dump(ref_object, open(ref_object.storage, "wb"))

	def load_reference_seqs(self):
		if len(os.listdir(self.ref_seq)) == 0:
			print(('no reference sequences stored in ' + self.ref_seq))

		else:

			for file in os.listdir(self.ref_seq):
				if '.pickle' in file:
					ref_path = os.path.join(self.ref_seq, file)
					reference = pickle.load(open(ref_path, "rb"))
					self.reference_seqs[file.replace(
						'.pickle', '')] = reference

	def store_reference_seqs(self):
		for key in self.reference_seqs:
			pickle.dump(self.reference_seqs[key], open(
				self.reference_seqs[key].storage, "wb"))
