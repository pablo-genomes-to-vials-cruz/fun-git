#modified from "http://enve-omics.ce.gatech.edu/enveomics/docs?t=FastA.N50.pl"
#now it also calculates GC content ;)
#now displays data in tabular format


open FILE, $ARGV[0] or die "i can not open the input FILE\n";
$name="$ARGV[0]";
while ($line=<FILE>){
    if ($line=~/>/){
    $contig_number++;    
    }
    
  else {
        $towaste="$line";
        $A+=$towaste=~s/A//gi;
        $C+=$towaste=~s/C//gi;
        $G+=$towaste=~s/G//gi;
        $T+=$towaste=~s/T//gi;
        $N+=$towaste=~s/N//gi;
        $length=$A+$T+$C+$G+$N;
       }
}

$GC_raw = eval { ($G+$C)/($A+$T+$C+$G)};
$GC=sprintf("%.2f",$GC_raw);
use List::Util qw/sum min max/;

my ($seqs, $minlen, $n__) = @ARGV;
$seqs or die "
Description:
   Calculates the N50 value of a set of sequences.  Alternatively, it
   can calculate other N** values.  It also calculates the total number
   of sequences and the total added length.
   
Usage:
   $0 seqs.fa[ minlen[ **]]
   seqs.fa	A FastA file containing the sequences.
   minlen	(optional) The minimum length to take into consideration.
   		By default: 0.
   **		Value N** to calculate.  By default: 50 (N50).
";
$minlen ||= 0;
$n__    ||= 50;

my @len = ();
open SEQ, "<", $seqs or die "Cannot open file: $seqs: $!\n";
while(<SEQ>){
   if(/^>/){
      push @len, 0;
   }else{


 next if /^;/;
      chomp;
      s/\W//g;
      $len[-1]+=length $_;
   }
}
close SEQ;
@len = sort { $a <=> $b } map { $_>=$minlen?$_:() } @len;
my $tot = (sum(@len) || 0);

my $thr = $n__*$tot/100;
my $pos = 0;
for(@len){
   $pos+= $_;
   if($pos>=$thr){
      $n50="$_";
      last;
   }
}
$GC_raw = eval { ($G+$C)/($A+$T+$C+$G)};
$GC=sprintf("%.2f",$GC_raw);
#print "\File_name\tAssembly_Length\tContig_number\tGC_content\tN50\n";
print "$name\t$tot\t$contig_number\t$GC\t$n50\t\n";
