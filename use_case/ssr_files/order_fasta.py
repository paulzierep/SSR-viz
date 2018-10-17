
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO

def get_keyword(rec):
	key = str(rec.id).split('_')[-1]
	return(key)


sorted_fasta = (sorted(list(SeqIO.parse("NRPS_A_dom_ali.fasta", "fasta")), key = get_keyword))

SeqIO.write(sorted_fasta, "NRPS_A_dom_ali_	sorted.fasta", 'fasta')

