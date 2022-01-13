#a simple SH wrapper to run augustus, then taking the gff file and using it to run antismash, then producting the tabular and fasta files for FUNGISON
#input is a .fna DNA fasta file, output is antimash file, gene calls *.gff protein fasta files, and tabular file
#name of fna file is the only input needed e. g. sh wrapper.sh Tolypocladium_inflatum.fna
LABEL=$(echo "$1" | sed 's/.fna//')
echo "running: augustus --species=Exophiala_dermatitidis  $LABEL.fna --gff3=on --stopCodonExcludedFromCDS=off > $LABEL.gff" 
/home/ynp/Augustus/bin/augustus --species=Exophiala_dermatitidis  $LABEL.fna --gff3=on --stopCodonExcludedFromCDS=off > $LABEL.gff
echo "running: antismash --taxon fungi --genefinding-gff3 $LABEL.gff $LABEL.fna" 
antismash --taxon fungi --fullhmmer --genefinding-gff3  $LABEL.gff $LABEL.fna
perl  GFF_GBK_to_FUNGISON_2.pl $LABEL.gff ./$LABEL/$LABEL.gbk
