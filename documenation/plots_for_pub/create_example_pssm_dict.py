from Bio.Seq import Seq
#from Bio.Alphabet import generic_protein
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

class_A = [	'AAAAALAAAAAIAAAAAAR',
			'AAAALLAAAAIIIIIIIAR',
			'AAALLLAAAIIILLLLLAR',
			'AALLLLAAIIIIMMMMMAR',
			'ALLLLLAIIIIIFFFFFAR',
		]

class_B = [	'LLLLLLLLLLLLLAAAALL',
			'LLLLLLLLLLLLLLQQQLL',
			'LLLLLLLLLLLLLLLWWLL',
			'LLLLLLLLLLLLLLLLCLL',
			'LLLLLLLLLLLLLLLLLLL',
]

discri_dict = {'A':[],'B':[]}

index = 0
for seq1, seq2 in zip(class_A, class_B):
	# sequence1 = Seq(seq1, generic_protein)
	# sequence2 = Seq(seq2, generic_protein)

	record1 = SeqRecord(Seq(seq1, IUPAC.protein), description=str(index))
	record2 = SeqRecord(Seq(seq2, IUPAC.protein), description=str(index))

	discri_dict['A'].append(record1)
	discri_dict['B'].append(record2)
	index += 1

print(discri_dict)


