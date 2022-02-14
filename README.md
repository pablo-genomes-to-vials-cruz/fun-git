# fun-git
A set of tools to take fungal genomes mine them for natural products, compare them using phylogenomics, cores, pangenomes. Fun As in Fungi

**Dependencies:**

AUGUSTUS (3.4.0)  a gene prediction tool.
Sources and documentation at https://github.com/Gaius-Augustus/Augustus

Antismash
Genome annotation tool specilaized in Natural products
https://docs.antismash.secondarymetabolites.org/install/

CORASON
The version in this repository
https://github.com/pablo-genomes-to-vials-cruz/myCORASON 

**How to use it:**
1. Sequence and assamble some genomes or download them. Use a fasta file

2. Define a taxonomic affiliation for your genome: E. g. Exophiala, create a gene model for it
Division: Ascomycota
Class: 	Chaetothyriomycetes
Order: 	Chaetothyriales
Family: 	Herpotrichiellaceae
Genus: 	Exophiala
Train a model for annotation (I have some premade already -see repository- ) detailed instructions are available in the Augustus manual
E.g the closest annotated reference genome in the Genbanks is  Exophiala dermatitidis

NW_008751646  Exophiala dermatitidis NIH/UT8656 
Therefore:
export AUGUSTUS_CONFIG_PATH=/home/ynp/Augustus/config 
perl new_species.pl --species=Exophiala dermatitidis
etraining --species=Exophiala_dermatitidis NW_008751646.gbk

3. Run Augustus to call the genes:
augustus --species=Exophiala_dermatitidis  mygenome.fna --gff3=on --stopCodonExcludedFromCDS=off > mygenome.gff
4. Run antismash to annotate the genes and find natural products BGCs
conda activate antismash (if you use conda to run antismash locally)
antismash --taxon fungi --fullhmmer --genefinding-gff3  mygenome.gff mygenome.fna


5. Integrate  the outputs into  CORASON-formatted files (the gbk file is the antismash whole genome file e.g. mysequence.gbk
perl  GFF_GBK_to_FUNGISON_2.pl mysequence.gff mysequence.gbk

6. Use this file to run corason

7. Enjoy



