from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import is_aa
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import os

# print(os.path.split('a/b/c.pdb'))

# exit()

def pdb2seq_2(pdbFile, chain, tmp):


	## First, open and parse the protein file
	p = PDBParser(PERMISSIVE=1)
	structure = p.get_structure(pdbFile, pdbFile)

	## Now go through the hierarchy of the PDB file
	##
	## 1- Structure
	##      2- Model
	##          3- Chains
	##              4- Residues
	##
	
	if (len(list(structure))) > 1:
		print('''Multiple structure parsing is not supported, plase create a pdb with only
			one structure''')
		exit()


	for model in structure:

		for pdb_chain in model:

			'''
			Seqres does start at the beginning of the protein, but often only a certain 
			part is shown in the pdb, therefore a certain number of - is added to the chain,
			so that the correct number corresponds to the sequence
			'''

			#get real starting number in the pdb
			number_in_pdb = ((pdb_chain.get_list()[0]).get_id()[1]) - 1

			#get gap to add to alignment
			head_attachment = ''.join(['-' for x in range(number_in_pdb)])

			seq = list()
			chainID = pdb_chain.get_id()

			if chain == chainID:

				for residue in pdb_chain:
					## The test below checks if the amino acid
					## is one of the 20 standard amino acids
					## Some proteins have "UNK" or "XXX", or other symbols
					## for missing or unknown residues
					if is_aa(residue.get_resname(), standard=True):
						seq.append(three_to_one(residue.get_resname()))
					else:
						seq.append("X")
				 
				## This line is used to display the sequence from each chain
				 
				print(">Chain_" + chainID + "\n" + str("".join(seq)))
				 
				## The next two lines create a sequence object/record
				## It might be useful if you want to do something with the sequence later
				
				final_seq = head_attachment + (str(''.join(seq)))

				myProt = Seq(final_seq, IUPAC.protein)
				record = SeqRecord(myProt, id=chainID, name="", description="")
				
				pdb_file_name = os.path.split(pdbFile)[-1].replace('.pdb', '')
				f_name = pdb_file_name + '_' + chain

				out_file = os.path.join(tmp, f_name + '.fasta')
				with open(out_file, "w") as output_handle:
					SeqIO.write(record, output_handle, "fasta")

				return(out_file)

		print('Chain not found in the {0} pdb file'.format(pdbFile))

## The end

# def pdb2seq(pdb, chain, tmp):
# 	'''
# 	Returns a sequence of the chain and pdb
# 	'''

# 	# print(chain)
# 	# print(pdb)
# 	# exit()

# 	handle = open(pdb, "rU")
	
# 	# print(handle)
# 	print(list(SeqIO.parse(handle, "pdb-seqres")))
# 	exit()

# 	for record in SeqIO.parse(handle, "pdb-seqres"):
# 		# print(record)
# 		# exit()
# 		chain_pdb = record.annotations['chain'] #get the chain

# 		if chain == chain_pdb:
# 			f_name = pdb + '_' + chain
# 			out_file = os.path.join(tmp, f_name + '.fasta')
# 			with open(out_file, "w") as output_handle:
# 				SeqIO.write(record, output_handle, "fasta")
# 			return(out_file)

# 	print('Chain not found in the {0} pdb file'.format(pdb))

#pdb = './pdbs/3czp_A.pdb'
pdb = './pdbs/5llb.pdb'

chain = 'A'
temp = './temp'

pdb2seq_2(pdb, chain, temp)
#pdb2seq(pdb, chain, temp)