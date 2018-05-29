import os
import subprocess
import pandas as pd

from Bio import SeqIO

from biopandas.pdb import PandasPdb

MAIN_FOLDER = '.'

#INPUT_PDB = '5ll0.pdb'
INPUT_PDB = '5lcd.pdb'
#INPUT_PDB = '1b3n.pdb'
#INPUT_PDB = '4na1.pdb'
#INPUT_ALI = 'phosphatases_ali.fasta'
INPUT_ALI = 'phosphatases_ali.fasta'

chain_id = 'A'

#positions = [650]	
#positions = [550,587,590,593,617,619,624,650,655,718]
# positions = [533,596]
# positions = [533,596] 

positions = [366,420] #I vs III
#positions = [459,496,533,536,560,563,565,570,596,664] #I vs III


#positions = [496,533,536,539,563,565,570,596,601,664]

pdb_file = os.path.join(MAIN_FOLDER, INPUT_PDB)
ali_file = os.path.join(MAIN_FOLDER, INPUT_ALI)

output_seq_file = INPUT_PDB.replace('.pdb', '_seq.fasta')
output_seq_file = os.path.join(MAIN_FOLDER, output_seq_file)

added_seq_alignment_file = INPUT_ALI.replace('.fasta', '') + '_' + INPUT_PDB.replace('.pdb','') + '.fasta'
added_seq_alignment_file = os.path.join(MAIN_FOLDER, added_seq_alignment_file)

mapoutput = output_seq_file + '.map' #thats how mafft does it !

pdb_name = INPUT_PDB.replace('.pdb','')

# Initialize a new PandasPdb object
# and fetch the PDB file from rcsb.org
ppdb = PandasPdb()
ppdb.read_pdb(pdb_file)
seq_df = ppdb.amino3to1()
seq_df['residue_number'] = ppdb.df['ATOM']['residue_number'] #map the sequence AS to the residue number


chain_seq = seq_df.loc[(seq_df.chain_id == chain_id)]
chain_seq_string = ''.join(list(chain_seq.residue_name))

pdb_seq_file = '>{0}\n{1}'.format(pdb_name, chain_seq_string)


if not os.path.exists(output_seq_file):

	with open(output_seq_file, 'w') as os_file:
		os_file.write(pdb_seq_file)

if not os.path.exists(mapoutput) or not os.path.exists(added_seq_alignment_file):

	stdout,stderr = subprocess.Popen(
		'mafft --anysymbol --addfull {0} --mapout --keeplength {1} > {2}'.format(
			output_seq_file,
			ali_file,
			added_seq_alignment_file),
			shell=True,
			stdout=subprocess.PIPE,
			stderr=subprocess.PIPE
			).communicate()

# print(stdout)
# print(stderr)

with open(mapoutput, 'r') as os_file:
	lines = list(os_file.readlines())

lines = lines[2:]
lines = [line.replace('\n','') for line in lines]
lines = [line.replace(' ','') for line in lines]
lines = [line.split(',') for line in lines]

labels = ['AS','pdb_seq', 'new_seq']

df_map = pd.DataFrame.from_records(lines, columns=labels)
chain_seq = chain_seq.reset_index()

mapped_alignment = pd.concat([df_map, chain_seq], axis = 1)
ali_res = mapped_alignment.loc[(mapped_alignment.new_seq != '-')]
ali_res = ali_res.residue_number


print('++++++++++++++++')
print(positions)
print('++++++++++++++++')

#########################
# Open with window
#########################

#!/usr/bin/env python
 
# Tell PyMOL to launch quiet (-q), fullscreen (-e) and without internal GUI (-i)
import __main__
__main__.pymol_argv = [ 'pymol', '' ]
 
import pymol
 
# Call the function below before using any PyMOL modules.
pymol.finish_launching()
 
from pymol import cmd

cmd.load(pdb_file)
cmd.hide('all')
cmd.show('cartoon', 'chain {0}'.format(chain_id))

# exit()
#resis = '+'.join(list(ali_res.astype('str')))

#resis = str(min(ali_res)) + '-' + str(max(ali_res))

selection_in_ali = positions
selection_in_ali = [str(num) for num in selection_in_ali]

# print('++++++++++++++++')
# print(mapped_alignment)
# print('++++++++++++++++')

selection = mapped_alignment.loc[(mapped_alignment.new_seq.isin(selection_in_ali))]

print('++++++++++++++++')
print(selection)
print('++++++++++++++++')

resis = '+'.join(list(selection.residue_number.astype('str')))

# ali_res = mapped_alignment.loc[(mapped_alignment.new_seq != '-')]


selection_str = 'chain {0} and resi {1}'.format(chain_id, resis)
print(selection_str)

cmd.select('aligned_area', selection_str)
cmd.color('red','aligned_area')
#cmd.show_as('cartoon', 'aligned_area')
#cmd.show_as('sticks', 'aligned_area')
cmd.zoom('aligned_area', buffer=15,state=0,complete=0)

cmd.show('cartoon', 'aligned_area')
cmd.show('sticks', 'aligned_area')
