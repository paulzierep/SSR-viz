from Bio.Seq import Seq
#from Bio.Alphabet import generic_protein
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

# class_A = [	'AAAAALAAAAAIAAAAAARDD',
# 			'AAAALLAAAAIIIIIIIARDD',
# 			'AAALLLAAAIIIGGGGGARDD',
# 			'AALLLLAAIIIIMMMMMARDD',
# 			'ALLLLLAIIIIIFFFFFARDD',
# 		]

# class_B = [	'LLLLLLLLLLLLLAAAALLNK',
# 			'LLLLLLLLLLLLLLQQQLLNK',
# 			'LLLLLLLLLLLLLLLWWLLNK',
# 			'LLLLLLLLLLLLLLLLCLLNK',
# 			'LLLLLLLLLLLLLLLLLLLNK',
# ]

class_A = [	'AAD',
			'AAD',
			'AAD',
			'AAD',
			'ALL']

class_B = [	'LLK',
			'LLK',
			'LLK',
			'LLK',
			'LLK',]




def get_example():
	discri_dict = {'A':[],'B':[]}

	index = 0
	for seq1, seq2 in zip(class_A, class_B):
		# sequence1 = Seq(seq1, generic_protein)
		# sequence2 = Seq(seq2, generic_protein)

		record1 = SeqRecord(Seq(seq1, IUPAC.protein), id = 'A', description=str(index))
		record2 = SeqRecord(Seq(seq2, IUPAC.protein), id = 'B', description=str(index))

		discri_dict['A'].append(record1)
		discri_dict['B'].append(record2)
		index += 1

	return(discri_dict)

# write to fasta

discri_dict = get_example()

all_recs = discri_dict['A'] + discri_dict['B']

from Bio import SeqIO

with open("Test_seqs.fasta", "w") as output_handle:
    SeqIO.write(all_recs , output_handle, "fasta")



# print(get_example())


