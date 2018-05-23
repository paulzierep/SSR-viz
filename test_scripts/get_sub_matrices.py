import pandas as pd
import numpy as np
from Bio.Data import IUPACData
import Bio.SubsMat.MatrixInfo as MI


def get_sub_matrix(name = 'basic', gap_importance = 1):

	#get an alphabet
	p_letters_gap = IUPACData.extended_protein_letters + '-'

	inverse_identity = 1 - pd.DataFrame(np.identity(len(p_letters_gap)),\
					index=list(p_letters_gap), \
					columns=list(p_letters_gap))

	if name == 'basic':
		#create basic matrix
		sub_matrix = pd.DataFrame(np.identity(len(p_letters_gap)),\
					index=list(p_letters_gap), \
					columns=list(p_letters_gap))

	else:
		#get an matrix
		p_letters_gap = IUPACData.extended_protein_letters + '-'
		empty_array = np.empty((len(p_letters_gap),len(p_letters_gap)))
		empty_array[:] = np.nan
		sub_matrix = pd.DataFrame(empty_array, \
						index=list(p_letters_gap), \
						columns=list(p_letters_gap))

		if name in MI.available_matrices:
			matrix = getattr(MI,name)
		else:
			print(('Matrix {0} does not exsist, chose one of those {1}'.format(name, MI.available_matrices)))

		# print(matrix)
		# exit()
		for AS, score in list(matrix.items()):

			if not pd.isnull(sub_matrix.loc[AS[0],AS[1]]) or not \
				pd.isnull(sub_matrix.loc[AS[1],AS[0]]): #simple check if matrix is not symetric
				print('Unsymmetric matrices are not supported at the moment')
				exit()

			if AS[0] == AS[1]:
				score = 0
				#print(AS[0])

			sub_matrix.loc[AS[0],AS[1]] = score
			sub_matrix.loc[AS[1],AS[0]] = score

	#normalization of entire dataframe
	sub_matrix = sub_matrix - sub_matrix.min().min()
	sub_matrix = sub_matrix / sub_matrix.max().max()

	#inverse dataframe
	sub_matrix = 1 - sub_matrix

	#remove the identity  A == A is useless
	sub_matrix = sub_matrix * inverse_identity

	#special ASs are not used at the moment
	sub_matrix.loc[['U','O','J'],:] = 0
	sub_matrix.loc[:,['U','O','J']] = 0

	#the importance of the gap can be set individually from the matrix
	sub_matrix.loc['-',:] = gap_importance
	sub_matrix.loc[:,'-'] = gap_importance
	sub_matrix.loc['-','-'] = 0

	#some matrices do not support all letters
	sub_matrix.fillna(value = 0, inplace = True)

	return(sub_matrix)

#matrix = get_sub_matrix(name = 'blosum62')

matrix = get_sub_matrix(name = 'pam30')
print(matrix.to_dict())

