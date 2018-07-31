#from Bio.Alphabet import IUPAC
import random
import pandas as pd
import pickle
import copy


from Bio.Data import IUPACData

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment

from .domain_extractor import clustalo_alignment_flat, align2discri_dict



'''some basic vars'''

p_letters = IUPACData.extended_protein_letters #the the letters which are used in biopyton.
p_letters_gap = IUPACData.extended_protein_letters + '-'

#p_letters = 'GCDAGCGC'

#p_letters 
#p_letters = p_letters_gap
# exit()

####################################
# Unit test functions
####################################

def get_random_seqs(num = 10, length = 100, difference = 0.6, initial_seq = None, names = None, random_seed = False):
	'''
	creates a biopyton sequence list of random sequences, very nice for debugging, difference sets the diff to the first random 
	created sequence, for the next seqs.
	'''

	if random_seed: 										#keep same random params for each run
		random.seed(1)

	if initial_seq:
		length = len(initial_seq)
	else:
		initial_seq =	''.join([p_letters[int(random.random() * len(p_letters))] for x in range(length)]) #totally random protein sequence

	seq_list = []
	for idx in range(num):
		new_seq = [x for x in initial_seq[:]]

		num_indices = int(round(difference * length)) #choice of how many position to alter
		if num_indices == 0:
			print('num of indices to alter is too small !!')
			break

		#print num_indices
		index_list = []
		counter = 0
		while True:
			index = int(random.random() * length) #get a random index
			if index in index_list:
				continue

			new_seq[index] = p_letters[int(random.random() * len(p_letters))] #change the random index

			index_list.append(index) #make sure each position is only altered ones
			counter += 1
			if counter == num_indices:
				break

		seq_sting = ''.join(new_seq)
		seq = Seq(seq_sting, IUPAC.extended_protein)

		if not names:
			names = 'test'

		record = SeqRecord(seq, id = 'Class_{0}_Idx_{1}'.format(str(names), str(idx)), name = 'Class_{0}_Idx_{1}'.format(str(names), str(idx)), description = 'for testing')
		seq_list.append(record)

	return(seq_list)


def get_random_discri_dict(num_classes = 2, num_seqs = 10, length = 100, seed_diff = 0.6, class_diff = 0.5, random_seed = None):

	'''
	creates a random discri dict for debugging,
	seed_diff --> 	initially the seed sequences are generated, which are used to
					create all the other classes, therefore this leeds to the diversity of the classes to each other
	class_diff --> 	Difference inside a class
	'''
	# print('num_classes', num_classes)
	# print('num_seqs', num_seqs)


	seed_seqs = get_random_seqs(num = num_classes, length = length, difference = seed_diff, random_seed = random_seed) #seeq sequences, should be quiet different

	#########################
	# uneven dist. of seqs per class
	#########################

	if type(num_seqs) == type([]): 
		if not len(num_seqs) == num_classes:
			print('num of seq args must match num of classes')
			exit()
		else:

			discri_list = []
			for discri_class, seq_num in zip(list(range(num_classes)), num_seqs):
				class_seqs = get_random_seqs(num = seq_num, length = length, names = str(discri_class), initial_seq = seed_seqs[discri_class], difference = class_diff)
				discri_list += class_seqs

	#########################
	# even dist. of seqs per class
	#########################

	else:

		discri_list = []
		for discri_class in range(num_classes):
			class_seqs = get_random_seqs(num = num_seqs, length = length, names = str(discri_class), initial_seq = seed_seqs[discri_class], difference = class_diff)
			discri_list += class_seqs

	alignment_path = clustalo_alignment_flat(INPUT_SEQ_LIST = discri_list, return_seqs = False, debug = False)
	discri_dict = align2discri_dict(alignment_path, regex_str = 'Class_([0-9]+)_')

	return(discri_dict)

####################################
# PSSM functions
####################################

def seq_list2pssm(seq_list, pandas = True, set_alphabet = True):
	'''
	creates a biopython/pandas pssm from a list of sequences (they should be aligned, logic for a pssm),
	'''

	seq_number = len(seq_list) #use the number of seqs to create the pwm
	align = MultipleSeqAlignment(seq_list)

	# the best way to avoid errors due to a wrong alphabet is to set it here, if the letter or gap does not exsist in the
	# alignment it does not really change anything, as this information would use drop in the mechine learning algorthms.
	# the biggest protein alphabet is used
	if set_alphabet:
		align._alphabet.letters = p_letters_gap


	summary_align = AlignInfo.SummaryInfo(align)
	consensus = summary_align.dumb_consensus()
	pssm = summary_align.pos_specific_score_matrix(consensus)

	# convert to pandas data frame
	if pandas:
		seq_list = [x for x in pssm] #get a list of dicts
		pssm = pd.DataFrame(seq_list)
		pwm = pssm.divide(seq_number) #position weight matrix (normalized pssm)
		return(pwm)

	# should be changed also to pwm or unconditional
	return(pssm)


def discri_dict2pssm_dict(discri_dict):

	'''
	convertes the entire discri dict, to a pssm dict, used for the pssm algo
	'''

	pssm_dict = {}
	for key, seqs in list(discri_dict.items()):
		pssm_dict[key] = seq_list2pssm(seqs)

	return(pssm_dict)

####################################
# Machine learning functions
####################################


def seq2OneHotEncoder(sequence):
	'''
	convertes a biopython sequence to a OneHotEncoded sequence
	'''

	letters = IUPAC.extended_protein.letters + '-'

	letter2digit = {}
	digit2letter = {}

	for num in range(0, len(letters)):											#create translation tables from letter to number
		letter2digit[sorted(letters)[num]] = num
		digit2letter[num] = sorted(letters)[num]

	seq_trans = [letter2digit[x] for x in sequence.seq] 						#translation as digits

	seq_OneHotEncoder = []														#create the binary list
	#label_names = []
	for num in seq_trans:
		binary_list = [0 if x != num else 1 for x in range(0, len(letters))]
		seq_OneHotEncoder += binary_list
		#label_name.append()

	return(seq_OneHotEncoder)

def get_labels_ML(sequence):

	'''
	Creates OneHotEncoded labels for a sequence
	'''

	letters = IUPAC.extended_protein.letters + '-'
	feature_names = []

	letter2digit = {}
	digit2letter = {}

	for num in range(0, len(letters)):											#create translation tables from letter to number
		letter2digit[sorted(letters)[num]] = num
		digit2letter[num] = sorted(letters)[num]

	for seq_pos in range(0, len(sequence.seq)):
		for letter in range(0, len(letters)):
			digit = digit2letter[letter] 
			feature_name = str(seq_pos) + '_' +  digit
			feature_names.append(feature_name)

	return(feature_names)


def discri_dict2ML_df(discri_dict):

	'''
	convertes the entire discri dict, to a pandas dataframe, 
	that stores the test_data and class label -> can then be used for all ML algos
	OneHotEncoded and Class annotated
	'''

	labels = get_labels_ML(discri_dict[list(discri_dict.keys())[0]][0]) 				#get the labels for the OneHotEncoded (from first seq in discri dict)
	labels += ['Class']															#store class ref in same dataframe

	seq_list = []
	for key, seqs in list(discri_dict.items()):
		for seq in seqs:
			seq_list.append(seq2OneHotEncoder(seq) + [key])

	df = pd.DataFrame(seq_list, columns = labels)								#create dataframe

	return(df)

'''
unit test for discri_dict2ML_df
'''

# discri_dict = get_random_discri_dict(num_classes = 4, num_seqs = 10, length = 100, seed_diff = 0.2, class_diff = 0.01)
# df = discri_dict2ML_df(discri_dict)

# print(discri_dict['1'][0].seq)
# df = df.iloc[0][:'0_J']
# print(df.to_frame().T.to_latex())

'''
unit test for the get_random_discri_dict function
'''

# discri_dict = get_random_discri_dict(num_classes = 4, num_seqs = 10, length = 100, seed_diff = 0.2, class_diff = 0.01)
# for keys, seqs in discri_dict.iteritems():
# 	print keys
# 	for seq in seqs:
# 		print(seq.seq)
# 		print(len(seq.seq))

'''
unit test for the get_random_seqs function
'''

# seqs = get_random_seqs(difference = 0.05, num = 10, length = 100)
# for seq in seqs:
# 	print seq.seq
#print seqs[0].seq

#print seq_list2pssm(seqs)

'''
unit test for the seq_list2pssm function
'''

# seqs = get_random_seqs(difference = 0.6)
# align = clustalo_alignment_flat(INPUT_SEQ_LIST = seqs, return_seqs = True, debug = False)# for seq in align:
# print seq_list2pssm(align)


'''
unit test for the discri_dict2pssm_dict function
'''

# discri_dict = get_random_discri_dict(num_classes = 2, num_seqs = 2, length = 30, seed_diff = 0.1, class_diff = 0.1)
# # for keys, seqs in discri_dict.iteritems():
# # 	print keys
# # 	for seq in seqs:
# # 		print(seq.seq)


# pssm_dict = discri_dict2pssm_dict(discri_dict)
# for keys, pssm in pssm_dict.iteritems():
# 	print keys
# 	print pssm





# def seq_list2onehotdecoder(seq_list):

# 	letter2digit = {}
# 	digit2letter = {}

# 	'''
# 	assigns a number to each letter in the alphabet
# 	'''

# 	for num in range(0, len(p_letters_gap)):
# 		letter2digit[sorted(p_letters_gap)[num]] = num
# 		digit2letter[num] = sorted(p_letters_gap)[num]

# 	data = []

# 	for seq in seq_list:
# 		seq_trans = [letter2digit[x] for x in seq.seq] #translation as digits
# 		seq_OneHotEncoder = []
# 		for num in seq_trans:
# 			binary_list = [0 if x != num else 1 for x in range(0, len(p_letters_gap))]
# 			seq_OneHotEncoder += binary_list

# 		data.append(seq_OneHotEncoder)

# 	return(data)

# random_discri_dict = get_random_discri_dict(num_classes = 3)










