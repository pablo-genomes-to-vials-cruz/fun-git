#creates a unique ID for this run
$rand=int(rand(1000000));
print "#########################################################\n";
print "## GBK GFF TO FUNGISON				       ##\n";
print "## A script to convert Augustus outputs in GFF3      ##\n";
print "## and antismash v6 (gbk) files into CORASON format    ##\n";
print "## by Pablo Cruz-Morales at DTU Biosustain 	       ##\n";
print "## pcruzm@biosustain.dtu.dk                	       ##\n";
print "## December 2021				       ##\n";
print "## Inputs: a GGF3 file with the AUGUSTUS gene calling   ##\n";
print "##        a GBK file with the  ANTISMASH 6 annotation  ##\n";
print "## usage: perl GFF_GBK_to_FUNGISON file.gff file.gbk    ##\n";
print "##run id : = $rand                                      ##\n";
print "#########################################################\n";
#Checking that the inputs are correct and providing instruction#
if ($ARGV[0]=~/.+.gff/){
	print "input Augustus gff3 file is $ARGV[0]\n";
	}
else {print "Error: did you entered a gff file?\nusage: perl GFF_GBK_to_FUNGISON file.gff file.gbk\n ";}

if ($ARGV[1]=~/.+.gbk/){
	print "input Antismash gbk file is $ARGV[1]\n";
	}
else {
	print "Error: did you entered a gbk file?\nusage: perl GFF_GBK_to_FUNGISON file.gff file.gbk\n ";
	}
#storing the file name for the outputs
open FILE, $ARGV[0] or die "I cant open the augustus GFF3 file\n";
$file_name_raw="$ARGV[0]";
$file_name_raw=~/(.+)(.gff)/;
$file_name="$1";
open TABLE,  ">$file_name.txt";
#making a temporary table file for protein sequences and ids RAST style with a random number for genome ID from the gff file
#watch out for different Augustus versions,  I used 3.4.0, differetn verions may have different comments, in such case update the following lines to get clean aminoacids
open AMINO,  ">amino_acids_$rand";
while ($line=<FILE>){
$line=~s/##gff-version 3//g;
$line=~s/# This output was generated with AUGUSTUS.+//g;
$line=~s/# AUGUSTUS is a gene prediction tool written.+//g;
$line=~s/# O. Keller, S. König, L. Gerischer, L. Romoth and Katharina Hoff.//g;
$line=~s/# Please cite.+//g;
$line=~s/# Using native and syntenically mapped cDNA alignments to improve de novo gene finding//g;
$line=~s/.+ioinformatics\/btn013//g;
$line=~s/.+insic information on sequences given.//g;
$line=~s/.+izing the parameters using config directory.+//g;
$line=~s/.+Using default transition matrix.+//g;
$line=~s/# admissible start codons.+//g;
$line=~s/.+ is in fasta format.//g;
$line=~s/.+the sequences in the input set.//g;
$line=~s/# ----- prediction on sequence number .+//g;
$line=~s/# Predicted genes for sequence number .+//g;
$line=~s/# \(none\)//g;
$line=~s/.+command line.+//g;
$line=~s/# augustus.+ //g;
$line=~s/# Evidence for.+//g;
$line=~s/# \% of transcript supported.+//g;
$line=~s/# CDS exons.+//g;
$line=~s/# CDS introns.+//g;
$line=~s/# 5'UTR exons and introns.+//g;
$line=~s/# 3'UTR exons and introns.+//g;
$line=~s/# hint groups fully obeyed.+//g;
$line=~s/# incompatible hint groups.+//g;
$line=~s/# Sources of extrinsic.+//g;
$line=~s/#     RM.+//g; 
$line=~s/.+\-\-gff3\=on \-\-stopCodonExcludedFromCDS\=off//g;
	if ($line=~/#/){
	chomp $line;
	$line=~s/\#//g;
	$line=~s/protein sequence \= \[/\t/g;
	$line=~s/start gene g/fig|666666.$rand.peg./g;
	$line=~s/.+end.+/\n/g;
	$line=~s/\]//g;
	$line=~s/ //g;
	print AMINO $line;
	}
}
close AMINO;
close FILE;
#making an amino acids fasta file from the amino_acids table
system "tr <amino_acids_$rand '\t' '\n'|sed 's/fig|/>fig|/' > $file_name.faa";
#Creating a functions table using a GBK file from antismash with fullhmmer option on 
open GEBEKA,  ">functions_table_$rand";
system  "grep \-E \'/gene\=\|\/description\' $ARGV[1]\| tr \'\\n\' \' \'\|sed \'s/\\/gene\/\\n/g\'\|sed \'s\/  \/\/g\' \> raw_grepped_$rand";
open FILEGREP, "raw_grepped_$rand";
while($line=<FILEGREP>){
	if ($line=~/\s+\n/){$dummy="1";next;}
	if ($line=~/"g\d+"\/description\=\".+/){
	$line=~/(="g\d+")(\/description\=\".+)/;
	$id="$1";
	$annotation="$2";
	$id=~s/\=\"//;
	$id=~s/"//;
	$annotation=~s/\/description\=\"//g;
	$annotation=~s/Condensation domain"/C-/g;
	$annotation=~s/Phosphopantetheine attachment site"/xCP-/g;
	$annotation=~s/AMP-binding enzyme C-terminal domain"AMP-binding enzyme"/A-/g;
	$annotation=~s/AMP-binding enzyme"/A-/g;
	$annotation=~s/KR domain"/KR-/g;
	$annotation=~s/Zinc-binding dehydrogenase"Alcohol dehydrogenase GroES-like domain"/ER-/g;
	$annotation=~s/Methyltransferase domain"/MT-/g;
	$annotation=~s/Polyketide synthase dehydratase"/DH-/;
	$annotation=~s/Acyl transferase domain"/AT-/g;
	$annotation=~s/Ketoacyl-synthetase C-terminal extension"Beta-ketoacyl synthase, C-terminal domain"Beta-ketoacyl synthase, N-terminal domain"/KS-/g;
	$annotation=~s/Beta-ketoacyl synthase, N-terminal domain"Beta-ketoacyl synthase, C-terminal domain"/KS-/g;
        $annotation=~s/N-terminal domain"Beta-ketoacyl synthase, C-terminal domain"Beta-ketoacyl synthase, N-terminal domain"/KS-/g;
	$annotation=~s/KS-Ketoacyl-synthetase C-terminal extension"/KS-/g;
	$annotation=~s/Zinc-binding dehydrogenase"/ER-/g;
	$annotation=~s/Male sterility protein"/Acyl-loading-/g;
	$annotation=~s/Alcohol dehydrogenase GroES-like domain"ER/ER-/g;
	$annotation=~s/Beta-ketoacyl synthase, C-terminal domain"Beta-ketoacyl synthase, N-terminal domain"/KS-/g;		
	$annotation=~s/Beta-ketoacyl synthase, N-terminal domain"/KS-/g;
	$annotation=~s/Starter unit:ACP transacylase in aflatoxin/Start_AT-/g;
	$annotation=~s/Thioesterase domain"/TE-/g;
	$annotation=~s/(.*)\1/$1/g;
	$annotation=~s/"/ /g;
	$annotation_line="$id\t$annotation\n";
	$annotation_line=~s/ \n/\n/g;
	$annotation_line=~s/\(RNA-dependent/RNA-dependent/g; 
	$annotation_line=~s/\(a.k.a. RRM, RBD, or/RRM/g;
	$annotation_line=~s/\// /g;
	$annotation_line=~s/\// /g;	
	$annotation_line=~s/\,/ /g;	
	print GEBEKA "$annotation_line";
	}
	else {
	chomp $line;	
	$id="$line";
	$annotation="unknown function";
	$id=~s/\=\"//;
	$id=~s/"//;
	print GEBEKA "$id\t$annotation\n"; 
	}
}
close FILEGREP;
close GEBEKA;
#getting a list of coordinates gene numbers and strand
system "grep transcript $ARGV[0]|grep '#' -v>transcript_list_$rand";
#making the table (txt) file
print TABLE "contig_id	feature_id	type	location	start	stop	strand	function	locus_tag	figfam	species	nucleotide_sequence	amino_acid	sequence_accession\n";
open COORDINATES, "transcript_list_$rand" or die "I cannot open the funtions table file\n";
while ($line=<COORDINATES>){
	if ($line=~/transcript/){
        $line=~/(.+)\t(AUGUSTUS)\t(transcript)\t(\d+)\t(\d+)\t(.+)\t(\+?\-?)\t(\.)\t(.+Parent=)(.+)/;
	$contig="$1";
	$up="$4";
	$down="$5";
	$strand="$7";
        $gene="$10";
        $gene=~s/g/fig|666666.$rand.peg./;
	$type="transcript";
	$location="chromosome";
	$locustag="$contig";
	$figfam="figfam";
	$species="$ARGV[0]";
	$nucleotide="ATCG";
	$accession="$locustag";
	#finding  the amino acid sequences fron the amino_acids file	
	open AMINOSEQ, "amino_acids_$rand" or die "i cant see the aminoacids table \n";
	while ($seq=<AMINOSEQ>){
		$seq=~/(.+)\t(.+)/;
		$label="$1";
		$aminoacid="$2";
		if ($label eq $gene){
	#finding  the annotations fron the functions_table file
	open FUNCTIONS, "functions_table_$rand" or die "i cant see the FUNCTIONS file \n";
		while ($id=<FUNCTIONS>){
			$id=~/(.+)\t(.+)/;
			$match="$1";
			$function="$2";
			$match=~s/g/fig|666666.$rand.peg./g;
			if ($match eq $gene){
	#checking the orientation of the gene, if negative then flipping the numbers, this makes arrrows go in the right direction
			if ($strand=~/\+/){
			$start="$up"; 
			$end="$down"; 
			}
			if ($strand=~/\-/){
			#else{
			$start="$down"; 
			$end="$up";
			}
			print  TABLE "$contig\t$gene\t$type\t$location\t$start\t$end\t$strand\t$function\t$locustag\t$figfam\t$species\t$nucleotide\t$aminoacid\t$accession\n";
			}
			else{
			$cont++;
			}
		}
		}
		else{
		$cont++;
		}
	}
	}
	else {
	$cont++;
	}
}
print TABLE "$contig\t$gene\t$type\t$location\t$start\t$end\t$strand\t$function\t$locustag\t$figfam\t$species\t$nucleotide\t$aminoacid\t$accession\n";
close COORDINATES;
close AMINOSEQ;
close COORDINATES;
close TABLE;
system "rm transcript_list_$rand raw_grepped_$rand amino_acids_$rand functions_table_$rand";
;;
