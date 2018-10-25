import csv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO

seqs = []

with open('aa_ID_code_all.csv') as csvfile:
     reader = csv.DictReader(csvfile, delimiter='\t')
     for row in reader:
         record = SeqRecord(Seq(row['Seq'], IUPAC.protein),
                    id=row['Name'] + '_spec_' +row['Spec'].lower(), name="",
                    description="")
         seqs.append(record)


SeqIO.write(seqs, "NRPS_A_dom.fasta", "fasta")


