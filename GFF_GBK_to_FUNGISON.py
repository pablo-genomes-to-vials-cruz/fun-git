import pandas as pd
import random
import sys
#creates a unique ID for this run
rand = random.randint(1,1000000)
print("#########################################################")
print("## GBK GFF TO FUNGISON				       ##")
print("## A script to convert Augustus outputs in GFF3        ##")
print("## and antismash v6 (gbk) files into CORASON format    ##")
print("## by Pablo Cruz-Morales at DTU Biosustain 	       ##")
print("## pcruzm@biosustain.dtu.dk                	       ##")
print("## December 2021				       ##")
print("## Inputs: a GGF3 file with the AUGUSTUS gene calling  ##")
print("##        a GBK file with the  ANTISMASH 6 annotation  ##")
print("## usage: python GFF_GBK_to_FUNGISON file.gff file.gbk ##")
print("##run id : = $rand                                     ##")
print("#########################################################")
#Checking that the inputs are correct and providing instruction#
GGF3_file = sys.argv[1] # File with AUGUSTUS gene calling
GBK_file = sys.argv[2] # File with antiSMASH 6 annotation
FileName = sys.argv[1][:-4]
# Storing the file name for the outputs
Ids,Sequences = [], []
contig_id_list, feature_id_list, start_list, stop_list, strand_list, amino_acid_list, sequence_accession_list, genes_list = [], [], [], [], [], [], [], []
with open(GGF3_file,'r',encoding = 'utf-8') as file:
	lines=file.readlines()
	for i,_ in enumerate(lines):
		if lines[i].startswith("# start gene"):
			iD = ">fig|666666." + str(rand) + ".peg." + str(lines[i].split(' ')[-1].split('g')[-1])
			Ids.append(iD.strip())
			sequence=''
			i+=1 # Increase line read
			contig_id_list.append(lines[i+1].split('\t')[0])
			feature_id_list.append(iD.strip()[1:])
			strand_list.append(lines[i+1].split('\t')[6])
			if lines[i+1].split('\t')[6] == '+':
				start_list.append(lines[i+1].split('\t')[3])
				stop_list.append(lines[i+1].split('\t')[4])
			else:
				start_list.append(lines[i+1].split('\t')[4])
				stop_list.append(lines[i+1].split('\t')[3])
		elif lines[i].startswith('# protein sequence'):
			while(lines[i].strip().startswith("# end gene") == False):
				sequence=sequence + lines[i].replace('#', '').strip()
				i+=1
			sequence = sequence.replace("protein", '').replace("sequence", '').replace("=", '').replace("[", '').replace("]", '').strip()
			Sequences.append(sequence.strip())
		else:
			i+=1
			pass
file.close()
# Making an amino acids fasta file from the amino_acids table
with open(FileName + '_test_file.faa', 'w') as file:
	for i in range(len(Sequences)):
		file.write(Ids[i] + '\n')
		file.write(Sequences[i] + '\n')
file.close()
#making the table (txt) file
function_dict={}
with open(GBK_file,'r',encoding = 'utf-8') as file:
	lines=file.readlines()
	for i,_ in enumerate(lines):
		if lines[i].startswith("                     /gene="):
			gene = lines[i].strip().split("=")[-1].split('g')[-1][:-1]
			genes_list.append(gene)
			function_dict[gene] = ''
		elif lines[i].startswith("                     /description="):
			FunctionFoundFlag = True
			function=lines[i].strip().split("=")[-1].replace('"', '')
			function_dict[gene] = function
		else:
			pass
file.close()
df = pd.DataFrame()
df['contig_id']=contig_id_list
df['feature_id']=feature_id_list
df['type']=["transcript"]*len(contig_id_list)
df['location']=["chromosome"]*len(contig_id_list)
df['start']=start_list
df['stop']=stop_list
df['strand']=strand_list
df['function']=function_dict.values()
df['locus_tag']=contig_id_list
df['figfam']=["figfam"]*len(contig_id_list)
df['species']=[GGF3_file]*len(contig_id_list)
df['nucleotide_sequence']=["ATCG"]*len(contig_id_list)
df['amino_acid']=Sequences
df['sequence_accession']=contig_id_list
df.to_csv(FileName + '_table_test_file.txt', sep="\t", index=False)