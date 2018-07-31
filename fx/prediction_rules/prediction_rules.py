'''the prediction rules should be designed in the following way:
input: desci_dict, output folder mapped to function, which can evaluate a new sequence of the same type
'''




#python native packages
import sys
import os
import subprocess
import pickle

#biopython
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio import SearchIO
from Bio.Alphabet import IUPAC

#important functions
from . import domain_extractor

#import sequence translation functions

'''
a lot of functions will be moved here from domain_extractor
'''

from . import sequence_translations

#basic functions
#sys.path.append('../')
from .basic import basic

#for the tree class
# from sklearn import tree
# from sklearn.ensemble import RandomForestClassifier

#import pydotplus 

#for the pssm class
#import pssm_functions

#for the pssm class new
from . import pssm_func


class rule_object():

	'''mother class for rule objects, init and setup are always the same'''

	def __init__(self, rule_name, rule_storage):

		'''
		set-up the needed paths for the object
		'''
		self.rule_name = rule_name
		self.rule_data_storage = os.path.join(rule_storage, self.rule_name)
		self.rule_storage = os.path.join(rule_storage, self.rule_name + '.rule')

	def setup_prediction_data(self, **kwargs):
		'''
		create folder and files to store the data needed to predict sequences
		'''
		if basic.check_dict_key_true(kwargs, 'keep_data_folder'):
			basic.create_folder(self.rule_data_storage, remove = False, new_folder = False)
		else:
			basic.create_folder(self.rule_data_storage, remove = True)

		basic.create_file(self.rule_storage, remove = True)

	def align_seq(self, seq2align):
		'''
		get a temporary alignment for the seq2predict 
		'''


		# print('\n*************')
		# print(seq2align)

		TEMP_ALIGNMENT_PATH = os.path.join(self.rule_data_storage , 'temp_alignment.fasta')
		self.test_seq_aligned = domain_extractor.add_seq2alignment(	INPUT_SINGLE_FASTA_PATH = seq2align,
															INPUT_FASTA_ALIGNED_PATH = self.alignment,
															OUTPUT_PATH = TEMP_ALIGNMENT_PATH, 
															temp_files = self.rule_data_storage,
															debug = False) 			#align input seq with multiple seq alignment,
																											#is important, so that positions are equal
		TEST_SEQ_ALIGNED_PATH = os.path.join(self.rule_data_storage , 'aligned_test_seq.fasta')


		# print(self.test_seq_aligned)
		# print(type(self.test_seq_aligned))
		# print('*************\n')

		SeqIO.write(self.test_seq_aligned, TEST_SEQ_ALIGNED_PATH, 'fasta')
		self.test_seq_aligned_path = TEST_SEQ_ALIGNED_PATH

	def store_alignment(self, discri_dict):
		'''
		To any prediction rule an alignment is needed to align the sequence one wants to predict.
		This is best stored with the rule
		'''
		ALIGNMENT_PATH = os.path.join(self.rule_data_storage , 'alignment.fasta')
		domain_extractor.discri_dict2alignment(discri_dict, ALIGNMENT_PATH) #as an alignemnt is needed for the prediction it is better stored as well
		self.alignment = ALIGNMENT_PATH

	def create_prediction_data(discri_dict, **kwargs):
		'''
		specific to each rule, stores the main information to predict a sequence
		'''
		pass

	def predict_sequence(self, seq_path, detail = False):
		'''
		specific to each rule, loads the infos needed for prediction and aligns the sequence to predict
		'''
		pass


# class ML_rule(rule_object):

# 	'''
# 	the tree class convertes a seqeunce alignment into the onehotencoder format, which can be used to create a decicion tree, 
# 	th format is needed as the sklearn algorithm can only work with true and false classification.
# 	onehotencoder example:
	
# 	ABCCC-->
# 	A=1, B=0, C=0 ; A=0, B=1, C=0 ; A=0, B=0, C=1 ; A=1, B=0, C=1 ; A=1, B=0, C=1
# 	'''

# 	def create_prediction_data(self, discri_dict, **kwargs):

# 		'''creates all the needed files for the prediction and initiates the prediction object, which is stored as pickle object'''

# 		graph, clf = dict2tree(discri_dict, **kwargs) #main tree function
# 		if graph:
# 			graph.write_pdf(os.path.join(self.rule_data_storage , 'entropy.pdf')) #graphical output

# 		self.clf = clf #store the classifyer object in this object

# 		if basic.check_dict_key_true(kwargs, 'plot_feature_importance'):
# 			print(list(clf.feature_importances_)) #todo get feature names !!!!

# 		self.store_alignment(discri_dict)

# 	def predict_sequence(self, seq_path, detail = False, **kwargs):

# 		#call the seq alignment method to create an alignment for the test sequence to the rule alignment
# 		self.align_seq(seq_path)

# 		OneHotEncoder = seq2OneHotEncoder(self.test_seq_aligned) #create a hotencoder for the seq

# 		best_match = self.clf.predict([OneHotEncoder])[0] #predict the sequence

# 		#print(best_match)
# 		probability = self.clf.predict_proba([OneHotEncoder])[0]
# 		#print(probability)
# 		classes_ordered = (self.clf.classes_) #classes as ordered in the tree object
# 		#print(classes_ordered)

# 		result_dict = {} #here more details for the tree could be implemented
# 		index = 0
# 		for num in probability:
# 			result_dict[classes_ordered[index]] = num
# 			index += 1

# 		if detail == True:
# 			return(best_match, result_dict)
# 		else:
# 			return(best_match)



# class tree_rule(rule_object):

# 	'''
# 	the tree class convertes a seqeunce alignment into the onehotencoder format, which can be used to create a decicion tree, 
# 	th format is needed as the sklearn algorithm can only work with true and false classification.
# 	onehotencoder example:
	
# 	ABCCC-->
# 	A=1, B=0, C=0 ; A=0, B=1, C=0 ; A=0, B=0, C=1 ; A=1, B=0, C=1 ; A=1, B=0, C=1
# 	'''

# 	def create_prediction_data(self, discri_dict, **kwargs):

# 		'''creates all the needed files for the prediction and initiates the prediction object, which is stored as pickle object'''

# 		graph, clf = dict2tree(discri_dict, **kwargs) #main tree function
# 		if graph:
# 			graph.write_pdf(os.path.join(self.rule_data_storage , 'entropy.pdf')) #graphical output

# 		self.clf = clf #store the classifyer object in this object

# 		if basic.check_dict_key_true(kwargs, 'plot_feature_importance'):
# 			print(list(clf.feature_importances_)) #todo get feature names !!!!

# 		self.store_alignment(discri_dict)

# 	def predict_sequence(self, seq_path, detail = False, **kwargs):

# 		#call the seq alignment method to create an alignment for the test sequence to the rule alignment
# 		self.align_seq(seq_path)

# 		OneHotEncoder = seq2OneHotEncoder(self.test_seq_aligned) #create a hotencoder for the seq

# 		best_match = self.clf.predict([OneHotEncoder])[0] #predict the sequence

# 		#print(best_match)
# 		probability = self.clf.predict_proba([OneHotEncoder])[0]
# 		#print(probability)
# 		classes_ordered = (self.clf.classes_) #classes as ordered in the tree object
# 		#print(classes_ordered)

# 		result_dict = {} #here more details for the tree could be implemented
# 		index = 0
# 		for num in probability:
# 			result_dict[classes_ordered[index]] = num
# 			index += 1

# 		if detail == True:
# 			return(best_match, result_dict)
# 		else:
# 			return(best_match)

# class pssm_rule(rule_object):

# 	'''this class creates a window of specific positions in an alignment, based on a specific sequence, the classification is done via 
# 	pssms for this window'''


# 	def create_prediction_data(self, discri_dict, **kwargs):
# 		'''creates all the needed files for the prediction and initiates the prediction object, which is stored as pickle object'''

# 		# if kwargs == {}:
# 		# 	#this mode uses the best pssms for the differens-pssms

# 		if basic.check_dict_key_true(kwargs, 'create_pssm_fasta_dd'):
# 		# if 'create_pssm_fasta_dd' in kwargs: #get a pssm as fasta for the discri dict
# 		# 	if kwargs['create_pssm_fasta_dd'] == True:
# 				pssm_functions.write_discri2pssm(discri_dict, self.rule_data_storage)

# 		pssm_dict = pssm_functions.descri_dict2pssm(discri_dict) #convert the alignments to pssms


# 		difference_table = pssm_functions.all_pssm(pssm_dict) 	# calculate the difference pandas table, ergo a table where the positions are shown, where
# 																# the one-to-one compared pssms differ the most.
# 		pssm_functions.draw_heatmap(difference_table, self.rule_data_storage)

# 		break_up = pssm_functions.specific_pssm_creator(difference_table, discri_dict, self.rule_data_storage, **kwargs)

# 		if break_up:
# 			print('PSSM could not be calculated, becouse there are not enought significant positions.')

# 		self.store_alignment(discri_dict)

# 	def predict_sequence(self, seq_path, detail = False):

# 		#get a temporary alignment
# 		TEMP_ALIGNMENT_PATH = os.path.join(self.rule_data_storage , 'temp_alignment.fasta')
# 		seq_a = domain_extractor.add_seq2alignment(seq_path, self.alignment, TEMP_ALIGNMENT_PATH) 	#align input seq with multiple seq alignment,
# 																									#is important, so that positions are equal

# 		result_dict = pssm_functions.score_pssm(seq_a, self.rule_data_storage)


# 		######################
# 		#todo: This does not cover the case where multiple best matches are possible! 
# 		######################

# 		max_score = max(result_dict.values()) #get the best match
# 		for match in result_dict:
# 			if result_dict[match] == max_score:
# 				best_match = match

# 		if detail == True:
# 			return(best_match, result_dict)
# 		else:
# 			return(best_match)


class pssm_rule_new(rule_object):

	'''this class creates a window of specific positions in an alignment, based on a specific sequence, the classification is done via 
	pssms for this window


	the new implementation uses numpy and pandas matrix frunctions for speed 
	and a updated function !
	'''


	def create_prediction_data(self, discri_dict, **kwargs):
		'''creates all the needed files for the prediction and initiates the prediction object, which is stored as pickle object'''

		pssm_dict =  sequence_translations.discri_dict2pssm_dict(discri_dict) 

		if basic.check_dict_key_not_none(kwargs, 'gap_imp'):
			gap_imp = kwargs['gap_imp']
		else:
			gap_imp = 1

		#this matrix is used for the final AS type evaluation, there the gap importance is 
		#always 0, as an exchange with a gap is little information at all
		if 'sub_matrix' in list(kwargs.keys()):
			sub_matrix = pssm_func.get_sub_matrix(name = kwargs['sub_matrix'], gap_importance = 0)
		else:
			sub_matrix = pssm_func.get_sub_matrix(name = 'basic', gap_importance = 0)

		#this is the initial matrix which is used to score the exchange of amino acids,
		#here the gap parameter can change the scoring significantly

		basic_sub_matrix = pssm_func.get_sub_matrix(name = 'basic', gap_importance = gap_imp)

		if 'm_weigth' in list(kwargs.keys()):
			matrix_weigth = kwargs['m_weigth']

		df = pssm_func.all_vs_all_pssm_cons(pssm_dict, 
											sub_matrix= sub_matrix, 
											b_sub_matrix = basic_sub_matrix, 
											matrix_weigth = matrix_weigth)	#get the cons_scores for all pssms
		self.df = df 															#add df to the class 
		self.pssm_dict = pssm_dict 												#add pssm dict to the class

		if basic.check_dict_key_true(kwargs, 'visual'):

			# optional plot name for pssm rule
			if 'plot_name' in kwargs:
				# if not '.pdf' in kwargs['plot_name']: #add .pdf if not specified
				# 	kwargs['plot_name'] += '.pdf'
				visual_path = os.path.join(self.rule_data_storage, kwargs['plot_name'])
			else:
				visual_path = os.path.join(self.rule_data_storage, 'temp_visual.pdf')

			df = pssm_func.add_all_vs_all(df) 										#compute the average for a nicer plot

			#Create plot stats
			# if basic.check_dict_key_true(kwargs, 'stats'):
			# 	pssm_func.get_stats(df = df, pssm_dict = pssm_dict, path = visual_path, **kwargs)

			# do not remove old plot by default
			#print(kwargs['delete_plot'])

			#tested
			if basic.check_dict_key_true(kwargs, 'delete_plot'):
				if basic.create_file(visual_path, remove = True):
					pssm_func.df2pssm_visual(df = df, path = visual_path, pssm_dict = pssm_dict, **kwargs)
			else:
				if not basic.create_file(visual_path, remove = False):
					print('Plot not created, chose different plot name, please !')
				else:
					pssm_func.df2pssm_visual(df = df, path = visual_path, pssm_dict = pssm_dict, **kwargs)
					
			# else:
			# 	pssm_func.df2pssm_visual(df = df, path = visual_path, **kwargs)

		# if basic.check_dict_key_true(kwargs, 'visual'):
		# 	visual_path = os.path.join(self.rule_data_storage, 'temp_visual.svg')
		# 	df = pssm_func.add_average(df) 										#compute the average for a nicer plot
		# 	pssm_func.df2pssm_visual(df, path = visual_path, **kwargs)

		#store the current alignment if needed
		if not basic.check_dict_key_true(kwargs, 'no_alignment'):
			self.store_alignment(discri_dict)

	def predict_sequence(self, seq_path, detail = False, **kwargs):

		if 'position_threshold' in kwargs:
			position_threshold = kwargs['position_threshold']
		else:
			position_threshold = 10

		#call the seq alignment method to create an alignment for the test sequence to the rule alignment
		self.align_seq(seq_path)
		scoring_dict = pssm_func.pssm_scoring(self.test_seq_aligned, self.pssm_dict)
		result_dict =  pssm_func.weighed_pssm_scoring(scoring_dict, weigth_matrix = self.df, no_weigth_matrix = False, position_threshold = position_threshold)

		max_score = max(result_dict.values()) #get the best match
		for match in result_dict:
			if result_dict[match] == max_score:
				best_match = match

		if detail == True:
			return(best_match, result_dict)
		else:
			return(best_match)



# class hmmer_rule(rule_object):
# 	'''class for the hmmer prediction rules, the design is intended, so that a new sequence is used as input and prediction is given as output
# 	the objects are stored as pickle object'''

# 	def create_prediction_data(self, discri_dict, **p_kwargs):

# 		'''creates all the needed files for the prediction and initiates the prediction object, which is stored as pickle object'''

# 		for discri in discri_dict:

# 			align_output = os.path.join(self.rule_data_storage , str(discri) + '.fasta')
# 			hmmer_output = os.path.join(self.rule_data_storage , str(discri) + '.hmm')
# 			#print hmmer_output

# 			align1 = MultipleSeqAlignment(discri_dict[discri])
# 			SeqIO.write(align1, align_output, "fasta")

# 			domain_extractor.hmm_build(align_output, hmmer_output)

# 		self.store_alignment(discri_dict)

# 	def predict_sequence(self, seq_path, detail = False, **kwargs):

# 		'''Performes a hmmer search for each of the hmmer files. Returns the best ranked result.'''

# 		if 'domT_threshold' in kwargs:
# 			domT_threshold = kwargs['domT_threshold']
# 		else:
# 			domT_threshold = 10


# 		self.align_seq(seq_path)

# 		result_dict = {}
# 		for file in os.listdir(self.rule_data_storage):

# 			if '.hmm' in file:

# 				hmm_file_path = os.path.join(self.rule_data_storage, file)
# 				hmm_search_path = os.path.join(self.rule_data_storage, file.replace('.hmm', '.search'))

# 				sub_report = subprocess.Popen('hmmsearch --domT={3} {0} {1} > {2}'.format(hmm_file_path, 
# 											self.test_seq_aligned_path,
# 											hmm_search_path,
# 											domT_threshold), 
# 											shell=True, 
# 											stdout=subprocess.PIPE, 
# 											stderr=subprocess.PIPE).communicate()
# 				if sub_report != ('',''):
# 					print(sub_report)

# 				record = SearchIO.read(hmm_search_path, "hmmer3-text")

# 				if len(record) > 1:
# 					print('Multiple record error !! for ' + file)
# 					exit()

# 				elif len(record) == 0:
# 					print('No record found, set score to 0')
# 					bitscore = 0
# 					result_dict[file.replace('.hmm', '')] = bitscore

# 				else:
# 					for rec in record:
# 						if len(rec) != 1:
# 							print('Multiple matches found, scores will be combined!!')

# 						bitscore_sum = [] 	#if there are multiple matches, this is due to a long gap (assuming that the seq. was preprocessed), therefore in this case
# 											# the scores are summed up
# 						for match in rec:
# 							bitscore_sum.append(match.bitscore)

# 						bitscore = sum(bitscore_sum)
# 						result_dict[file.replace('.hmm', '')] = bitscore

# 		#print result_dict


# 		max_score = max(result_dict.values()) #get the best match
# 		for match in result_dict:
# 			if result_dict[match] == max_score:
# 				best_match = match

# 		if detail == True:
# 			return(best_match, result_dict)
# 		else:
# 			return(best_match)


def seq2OneHotEncoder(sequence):

	'''todo'''

	letters = IUPAC.extended_protein.letters + '-'

	letter2digit = {}
	digit2letter = {}

	for num in range(0, len(letters)):
		letter2digit[sorted(letters)[num]] = num
		digit2letter[num] = sorted(letters)[num]

	seq_trans = [letter2digit[x] for x in sequence.seq] #translation as digits
	seq_OneHotEncoder = []

	for num in seq_trans:
		binary_list = [0 if x != num else 1 for x in range(0, len(letters))]
		seq_OneHotEncoder += binary_list

	return(seq_OneHotEncoder)

def dict2tree(discri_dict, **kwargs):

	##################################################
	# create dict to translate AS code to digits
	##################################################

	# print(discri_dict[discri_dict.keys()[0]][0].seq.alphabet.letters)
	# exit()

	letters = IUPAC.extended_protein.letters + '-'

	seq_lenght = len(discri_dict[list(discri_dict.keys())[0]][0].seq)

	letter2digit = {}
	digit2letter = {}

	for num in range(0, len(letters)):
		letter2digit[sorted(letters)[num]] = num
		digit2letter[num] = sorted(letters)[num]
		

	##################################################
	#create data table in sklearn format
	##################################################

	'''sklearn cannot work with discret numbers yet, but no problem,
	the data just needs to be designed in a workable manner: Create one list for each position in the seq.
	then fill the list with binarys for each As code, the final result will be a tree where the 
	discriptor is the probability of a specific AS at a certain position.
	See here:
	https://datascience.stackexchange.com/questions/5226/strings-as-features-in-decision-tree-random-forest
	'''

	data = [] #data used to create the tree, must be a digit
	target = [] #descriptor

	#print discri_dict

	for descri in discri_dict:
		for seq in discri_dict[descri]:
			seq_trans = [letter2digit[x] for x in seq.seq] #translation as digits
			seq_OneHotEncoder = []
			for num in seq_trans:
				binary_list = [0 if x != num else 1 for x in range(0, len(letters))]
				seq_OneHotEncoder += binary_list

			#print seq_OneHotEncoder

			#if descri == 20:
			# 	print seq
			# 	print seq_OneHotEncoder

			data.append(seq_OneHotEncoder)
			target.append(descri)

	#print data
	#exit()
	# print target

	##################################################
	#define the class names and features
	##################################################

	feature_names = []

	for num1 in range(0, seq_lenght):
		for num2 in range(0, len(letters)):
			l = digit2letter[num2] 
			feature_name = str(num1) + l
			feature_names.append(feature_name)

	class_names = [str(x) for x in list(discri_dict.keys())]


	##################################################
	#build tree using sklearn
	##################################################

	if basic.check_dict_key_true(kwargs, 'RandomForrest'): #use the random forrest
			if 'algo_dict' in kwargs:
				clf = RandomForestClassifier(**kwargs['algo_dict']) #pass on the arguments to the function 
			else:
				clf = RandomForestClassifier()

			clf = clf.fit(data, target)

			#print(data)

			#for a,b in zip(data.columns, clf.feature_importances_):
				#print(a, b)


			graph = None #random forrest do not have one single tree
			return(graph, clf)

	else:
		#print(kwargs)
		if 'algo_dict' in kwargs:
			#print(kwargs['algo_dict'])
			clf = tree.DecisionTreeClassifier(**kwargs['algo_dict'])  #pass on the arguments to the function 
		else:
			clf = tree.DecisionTreeClassifier()

		clf = clf.fit(data, target)

		if basic.check_dict_key_true(kwargs, 'tree_fig'):

			dot_data = tree.export_graphviz(clf, 
			out_file=None,
			node_ids=True,
			class_names= class_names,
			feature_names = feature_names,
			)

			graph = pydotplus.graph_from_dot_data(dot_data) 

		else: 
			graph = None

		return(graph, clf)





def prediction_dispatcher(name):
	'''returns the function corresponding to the name, new prediction functions must be registered here'''
	prediction_dict = {
						#'hmmer': hmmer_rule,
						#'tree': tree_rule,
						#'pssm': pssm_rule,
						'pssm_new': pssm_rule_new,
						'test': lambda discri_dict, output: print(discri_dict, output),
						}

	if name in prediction_dict:
		return(prediction_dict[name])
	else:
		print(name + ' no such function is registered in the prediction_dispatcher')

