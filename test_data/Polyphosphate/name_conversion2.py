from Bio import SeqIO

records1 = list(SeqIO.parse('./phosphatases_ali_5ll0.fasta', "fasta"))

#records2 = list(SeqIO.parse('./promals3d_alignment.fasta', "fasta"))

records1 = sorted(records1, key = lambda seq: seq.description[-3:-1])
SeqIO.write(records1, './phosphatases_ali_5ll0_sorted.fasta', 'fasta')
#print(records1)

# for rec1 in records1:
# 	print(rec1)

	#print(rec1.id)
# 	for rec2 in records2:
# 		#print(rec2.description)
# 		#exit()
# 		if rec1.id in rec2.description:
# 			rec2.id = rec1.description
# 			rec2.name = ''
# 			rec2.description = ''

# for seq in records2:
# 	print(seq)

# SeqIO.write(records2, './promals3d_alignment_trim.fasta', 'fasta')