# For more information or help visit http://bioinfo.ut.ee/strainseeker

use strict;
use warnings;
use Getopt::Long;

################################
 ####### TIME MANAGEMENT ######
################################
my $prevTime = time;
my @timeArray;

################################
 ############ INPUT ###########
################################
# Get StrainSeeker path
my $path;
if (index($0, "/") == -1){
	$path = "./";
} else {
	my @pathArray = split("/",$0);
	pop(@pathArray);
	$path = (join ("/", @pathArray))."/";
}

# Options stuff
my $tree_data_final; # Database name
my @sample_fasta;
my $help;
my $projectname = "StrainSeeker_output";
my $projectnameclean = $projectname;
my $suffix = localtime($^T);
$suffix =~ s/[\s:]+//g;
my $seekerTemp = "strainseekertemp".$suffix;
my $gc_seekerTemp = "gc_seekerTemp".$suffix;
# Constants
my $s_lim = 0.025; #min % in node or strain to consider
my $pval_limit = 0.05/5000; # P-value cutoff divided by average number of comparisons. If OE pval is larger, then OE = 1.
my $min_total_kmers = 50; # Minimum number of k-mers in nodes. If subtree is being analyzed, jump over these nodes
my $version;
my $verbose = 0;
GetOptions(
    'd=s' => \$tree_data_final,
    'dir=s' => \$tree_data_final,
    'i=s@' => \@sample_fasta,
    'h' => \$help,
    'help' => \$help,
    'v' => \$version,
    'version' => \$version,
    'o=s' => \$projectname,
    'output=s' => \$projectname,
    #'m=f' => \$s_lim,
    #'min=f' => \$s_lim,
    #'p=f' => \$pval_limit,
    #'pval=f' => \$pval_limit,
    #'k=i' => \$min_total_kmers,
    #'mintot=i' => \$min_total_kmers,
    'verbose' => \$verbose
     ) or die printHelp()."\n";

if (!@sample_fasta){
	@sample_fasta = @ARGV;
}

if ($version){
	die "Seeker v1.5\n";
}

if ((!$tree_data_final) || (!@sample_fasta) || ($help)){
	die printHelp()."\n";
}

if (substr($tree_data_final, -1) eq "/"){
	chop $tree_data_final;
}

my $dir_tree = $tree_data_final;
my $DbBinary = $dir_tree."/db_binary";
$tree_data_final .= "/info.txt";
my %treeinf; #Tree info file hash (each key has array: [self][parent][child1][child2][count])
my @smalltrees; #Subroots from info file (array with names only)
my $word; #Word size for tree database, from info file
my $n_rootleaf = 0;

foreach (@sample_fasta) {unless (-e $_) { die("File $_ does not exist!\n"); }}	# Check if fastas exist
unless (-e $tree_data_final) { die("File $tree_data_final doesn't exist!\n"); } # Check if info file exist
unless (-d $dir_tree) { die("Directory $dir_tree doesn't exist!\n"); } # Check if dir with list files exist

readInfo();

################################
 ##### GLOBAL PARAMETERS ######
################################
# Coverage and other constants
my $cov_blim = 0.1;
my $cov_ulim = 0.9;

my $subroot_constant = 1;
my $current_subroot;

# GenomeTester4 programs and R-scripts paths
my $glm = $path."glistmaker"; # GListmaker
my $glc = $path."glistcompare"; # Glistcompare
my $glq = $path."glistquery"; # Glistquery
my $gc  = $path."gmer_counter"; #GmerCounter
my $oe_script = $path."oe.R";
my $cov_script = $path."cov.R";

# File locations
#my $whitelist = "$dir_tree/whitelist_$word\.list"; #Union of all tree k-mer lists (DB)
#my $subWhite = "$dir_tree/subwhite_$word\.list"; #Union of all subtree root lists (DB)
my $sample = "$seekerTemp\_$word\.list"; #Sample already intersected with whitelist (TMP)
my $tmp_subroots = "strainseeker_subroots".$^T; #Temp dir for keeping subroot counts (TMP)
my $outputfilename = $projectname;

my %samplehas; # KEY: node, VALUE: found in that node, hash used by sub inSample
my %sampleDistribution;
my $outputcount = 0; # Number of strains found in sample
my %finalout; # Final output hash
my @warnings; # Here are warnings stored
my $readlen = getReadLen($sample_fasta[0]);

# Strings
my $known = "KNOWN";
my $unknown = "RELATED";

### END GLOBALS ###


################################
 ####### FILE MANAGEMENT ######
################################

my @fileMustNotExist = ("sampletemp.txt", "nodetemp.txt", $sample,
						"subwhite\_$sample, temp_$word\_intrsec.list",
						"$seekerTemp\_$word\_intersec\.list");
my @dirMustNotExist = ("tmp_subroots");

foreach (@fileMustNotExist){
	if (-e $_) {die "Please remove $_ from working dir\n"}
}
foreach (@dirMustNotExist){
	if (-d $_) {die "Please remove $_ from working dir\n"}
}

################################
 ############# SUBS ###########
################################

sub add2time {
	my $out = "$_[0]\t".(time - $prevTime);
	push (@timeArray, $out);
	$prevTime = time;
}

sub covered {
	(my $current) = @_;

	my $total = 0;

	my $lines = 0;
	foreach (@{$sampleDistribution{$current}}){
		$total += $_;
		last if ($_ == 0 && $total != 0);
		$lines++;
	}
	if ($total == 0){
		push (@warnings, "WARNING: total in covered sub is 0, returning 1, coverages are not trustworthy");
		return 1;
	}

	my $total2 = $treeinf{$current}[4];

	if ($total2 == 0){
		push (@warnings, "WARNING: total2 in covered sub is 0, returning 1, coverages are not trustworthy");
		return 1;
	}

	# SAMPLE
	my $b_point;
	my $b_point_counter = 0;
	my $i = 1;
	foreach (@{$sampleDistribution{$current}}){
		$b_point_counter+=($_/$total);
		if ($b_point_counter >= $cov_blim){
			$b_point = $i;
			last;
		}
		$i++;
	}

	$i = $total;
	my $u_point;
	my $u_point_counter = 0;
	my $switch = 1;
	foreach (reverse @{$sampleDistribution{$current}}){
		next if $_ == 0;
		$u_point_counter+=($_/$total);
		$u_point = $lines;
		$lines--;
		if ($u_point_counter <= (1-$cov_ulim) || $lines == 0){
			last;
		}
		last if ($switch && $u_point_counter > (1-$cov_ulim));
		my $switch = 0;
	}
	$i = 0;
	my $good_total = 0;
	my $sum = 0;
	foreach (@{$sampleDistribution{$current}}){
		$i++;
		if ($b_point > $i){
			next;
		} elsif ($u_point < $i){
			last;
		}
		$sum+=$i*$_;
		$good_total+=$_;
	}

	if ($good_total == 0){
		push (@warnings, "WARNING: good_total in covered sub is 0, returning 1, coverages are not trustworthy");
		return 1;
	}

	my $mean = $sum/$good_total;

	my @calc_cov = qx/Rscript $cov_script $mean $b_point $u_point 1 $total2 2> \/dev\/null/;
	return (split(" ", $calc_cov[0]))[1];
}

sub checkSubroots { # Check all subroots and return the ones that have enough kmers. Argument - array of all subroots. Return - array of subroots above cutoff
	my @returnarray;
	my @allsubroots = @_;

	foreach (@allsubroots) {
		if (inSample($_)/$treeinf{$_}[4] > $s_lim) {
			push(@returnarray, $_);
		}
	}
	return @returnarray;
}

sub dynamicLim { # For dynamic min k-mer count calculation. Input - node name and k-mer amount in subroot. Return - minimum percentage of k-mers
	my $count = $treeinf{$_[0]}[4];
	my $constant = $_[1];
	#return $s_lim;
	return ((73/sqrt($count)+$constant)/100);
}

sub findOE { # Finds whether O/E is 1 (returns 0) or lesser(-1)/greater(1). # argument - node name
	my ($node) = @_;
	my $child1 = $treeinf{$node}[2];
	my $child2 = $treeinf{$node}[3];
	my ($n_liik1, $n_liik2, $n_koos, $l_liik1, $l_liik2, $l_koos) = ($treeinf{$child1}[4], $treeinf{$child2}[4], $treeinf{$node}[4], inSample($child1), inSample($child2), inSample($node));
	my @out = qx/Rscript $oe_script $n_liik1 $n_liik2 $n_koos $l_liik1 $l_liik2 $l_koos $readlen $word/;
	my $oe = (split(/\s+/, $out[1]))[1];
	my $pval = (split(/\s+/, $out[4]))[1];
	print STDERR "--> OE = $oe | p-val = $pval | C1 found: ".inSample($child1)."/$treeinf{$child1}[4] | C2 found: ".inSample($child2)."/$treeinf{$child2}[4]\n" if $verbose;
	if ($pval > $pval_limit) {
		print STDERR "Continue along this branch...\n" if $verbose;
		return 0;
	} else {
		if (($oe > 0.85) && ($oe < 50)) {
			print STDERR "Print out closest relatives\n" if $verbose;
			return 1;
		}
		print STDERR "($treeinf{$child1}[4], $treeinf{$child2}[4], $treeinf{$node}[4], inSample($child1), inSample($child2), inSample($node)\n" if $verbose;
		print STDERR "Node: $node has EXCEPTIONALLY low/high OE $oe\n" if $verbose;
		return -1;
	}
}

sub getKmerCounts{ # Populate %samplehas
	my $gc_in = join(" ", @sample_fasta);
	system ("$gc -dbb $dir_tree/db_binary --unique --distribution 1000 $gc_in > $gc_seekerTemp");

	open(GCFILE, $gc_seekerTemp) or die "Could not open $gc_seekerTemp";

	while(<GCFILE>){
		chomp;
		my @arr = split(/\t/, $_);
		my $name = shift @arr;
		shift @arr;
		$samplehas{$name} = shift @arr;
		shift @arr;
		@{$sampleDistribution{$name}} = @arr;
	}

	close GCFILE;
}

sub getReadLen {
	open (FASTA,$_[0]) or die "Cannot open fasta file for getting read length\n";
	my $i = 1;
	my $typeChar;
	my $lenString;
	while(<FASTA>){
		chomp;
		if (!$typeChar){
			$typeChar = substr($_, 0, 1);
			next;
		}

		if ($typeChar eq substr($_, 0, 1) || substr($_, 0, 1) eq "+"){
			last;
		}
		$lenString .= $_;

		$i++;
	}

	return length($lenString);
}

sub inSample { # Argument as node, return found in node count
  (my $node) = @_;
  if ($treeinf{$node}[4] == 0){
  	return 0;
  }
  if (exists $samplehas{$node}) {
      return $samplehas{$node};
  } else {
	  return 0;
	}
}

sub lastOKnode {
	(my $current) = @_;
	my $current_found = inSample($current);
	my $parent = $treeinf{$current}[1];
	if ($treeinf{$current}[4] < $min_total_kmers){ # If k-mer count is too little then search further
		lastOKnode ($parent);
	} else {
		my $current_found_ratio = $current_found/$treeinf{$current}[4]; # %found in current
		if ($current_found_ratio > $s_lim){
			return $current;
		} else {
			lastOKnode ($parent);
		}
	}
}

sub printChildren { # Print all children for given node
	(my $current) = @_;
	my $child1 = $treeinf{$current}[2];
	my $child2 = $treeinf{$current}[3];

	my $c1_percent = 0; # check if one child's found % is twice of other child's %, true => that, false => both
	my $c2_percent = 0;
	if ($treeinf{$child1}[4]){
		$c1_percent = inSample($child1)/$treeinf{$child1}[4];
	}
	if ($treeinf{$child2}[4]){
		$c2_percent = inSample($child2)/$treeinf{$child2}[4];
	}

	if ($c1_percent/2 > $c2_percent){
		return printChildren2($child1);
	}
	elsif ($c2_percent/2 > $c1_percent){
		return printChildren2($child2);
	}
	else {
		return printChildren2($child1).printChildren2($child2);
	}
}

sub printChildren2 {
	(my $current) = @_;
	my $child1 = $treeinf{$current}[2];
	my $child2 = $treeinf{$current}[3];

	if ($child1 eq "NA"){ # if current is strain
		my $found_ratio = 0;
		if ($treeinf{$current}[4] > 0){
			$found_ratio = inSample($current)/$treeinf{$current}[4];
		}
		return "$current\t".$found_ratio."\n";
	}
	else {
		return printChildren2($child1).printChildren2($child2);
	}
}

sub printHelp {
	print "Usage: $0 -d <DB DIR NAME> -i <SAMPLE.fastq> [OPTIONAL PARAMETERS]\n";
	print "Options:
	-h, --help\t - Print this help
	-v, --version\t - Print version of the program
	-i, none\t - Input file (can be multiple, each with own flag)
	-o, --output\t - Output file name (default $projectname)
	-d, --dir\t - Path to database directory
	-verbose\t - Print out more of the working process\n";
	return "";
}

sub printStrains {
	(my $current, my $lastOK) = @_;
	my $child1 = $treeinf{$current}[2];
	my $child2 = $treeinf{$current}[3];

	my $current_found = inSample($current); # number of node/strain k-mers found in sample
	if ($current_found){
		my $current_found_ratio = $current_found/$treeinf{$current}[4]; # %found in current
		print STDERR ("\nNode: $current lastOK: $lastOK current found: $current_found current total: $treeinf{$current}[4] found %: ") if ($verbose);
		print STDERR sprintf("%.2f\n",$current_found_ratio*100) if ($verbose);
	} else {
		print STDERR ("\nNode: $current lastOK: $lastOK current found: 0 current total: $treeinf{$current}[4] found %: 0\n") if ($verbose);
	}

	if (($child1 eq "NA") && ($child2 eq "NA")) { # Is strain
		if ($treeinf{$current}[4] < $min_total_kmers) {
			return 0;
		} else {
			return (checkModule($current, 1));
		}
	} else {
		return skipNode($current, $lastOK); # Is node
	}
	return 1;
}

sub checkModule {
	(my $current, my $type) = @_; # Type = strain(1) or node(0)
	my $current_found = inSample($current); # number of node/strain k-mers found in sample
	my $current_found_ratio = $current_found/$treeinf{$current}[4]; # %found in current

	if ($current_found_ratio > dynamicLim($current, $subroot_constant)) { # HAS TO BE MORE THAN THAT TO EVEN CONSIDER
		if ($type) { # Is strain
			#print self
			@{$finalout{$current}} = ($known,$current,covered($current));
			print STDERR "# KNOWN $current\n-----\n" if $verbose;
			return 1;
		} else { # Is node
			my $oe = findOE($current);
			if ($oe == 0) { # IF IS BETWEEN LIMITS (STRAIN UNDER THIS NODE)
				my $child1_return = printStrains($treeinf{$current}[2], $current);
				my $child2_return = printStrains($treeinf{$current}[3], $current);
				my $subroot_ratio = inSample($current_subroot)/$treeinf{$current_subroot}[4];
				if (!$child1_return && !$child2_return && ($subroot_ratio > dynamicLim($current_subroot, 7))) {
					# Print current
					@{$finalout{$current}} = ($unknown,printChildren($current),covered($current));
					print STDERR "# RELATED TO GROUP" if $verbose;
					foreach ($finalout{$current}[1]){
						print STDERR "\n$_" if $verbose;
					}
					print STDERR "\n-----\n" if $verbose;
				}
				return 1;
			}	elsif ($oe == 1) { # IF IS HIGHER THAN UPPER LIMIT: NEW STRAIN
					@{$finalout{$current}} = ($unknown,printChildren($current),covered($current));
					print STDERR "# RELATED TO GROUP" if $verbose;
					foreach ($finalout{$current}[1]){
						print STDERR "\n$_" if $verbose;
					}
					print STDERR "\n-----\n" if $verbose;
					return 1;
			} else {
				return 0;
			}
		}
	}
	return 0;
}

sub skipNode {
	(my $current, my $lastOK) = @_;
	my $child1 = $treeinf{$current}[2];
	my $child2 = $treeinf{$current}[3];
	if ($treeinf{$current}[4] < $min_total_kmers) { # Skip this node
		my $child1_return = printStrains($child1, $lastOK);
		my $child2_return = printStrains($child2, $lastOK);
		my $sum = $child1_return + $child2_return;
		return $sum;
	} else { # Check this node
		return checkModule($current, 0);
	}
}

sub readInfo {
	# Create hash, current node as key, whole line in array as value, after "---" get subtree roots
	print STDERR "Read tree info file...\n";
	my $isSubtree = 0; #Switch for reading subroots

	open (TREE, '<', $tree_data_final) or die("Tree info file $tree_data_final is missing!");
	while (<TREE>){
		chomp;
		if ($_ eq "---"){
			$isSubtree++;
			next;
		}
		if (!$isSubtree) { # Read info file header
			my @line = split("\t",$_);
			$word = (split(":",$line[1]))[1];
		} elsif ($isSubtree == 1) { # Read leaves + nodes
			my @line = split("\t",$_);
			chomp $line[4];
			@{$treeinf{$line[0]}} = @line;
			#if(!(-e "$dir_tree/$line[0]\_$word\.list")) {
			#	die("FILE $dir_tree/$line[0]\_$word\.list is MISSING!!!");
			#}
		} else { # Read subroots
			if ($_ ne "") {push(@smalltrees, $_);}
			if ($treeinf{$_}[2] eq "NA"){$n_rootleaf++;}
			next;
		}
	}
	close TREE;
	add2time("Reading tree info file");
}

sub subrootConstant { #If subroot coverage is very high, modify the minimum allowed percentage of k-mers found. Input - name of subroot. Returns the constant.
	my $count = $treeinf{$_[0]}[4];
	my $current_found = inSample($_[0]);
	my $percent = $current_found/$count;
	if ($percent < 0.5) {
		return 1;
	} elsif ($percent < 0.6) {
		return 2;
	} elsif ($percent < 0.7) {
		return 3;
	} elsif ($percent < 0.8) {
		return 4;
	} elsif ($percent < 0.9) {
		return 5;
	} elsif ($percent < 0.96) {
		return 7;
	} else {
		return 10;
	}
}

sub sec2readable {
	my $secs = shift;
	if    ($secs >= 365*24*60*60) { return sprintf '%.1fy', $secs/(365*24*60*60) }
	elsif ($secs >=     24*60*60) { return sprintf '%.1fd', $secs/(24*60*60) }
	elsif ($secs >=        60*60) { return sprintf '%.1fh', $secs/(60*60) }
	elsif ($secs >=           60) { return sprintf '%.1fm', $secs/(60) }
	else                          { return sprintf '%.1fs', $secs}
}

################################
 ############# MAIN ###########
################################
#prepSample(); # Conver sample to k-mers, compare with whitelist and subwhite
print STDERR "Getting k-mer counts from sample...\n";
getKmerCounts();
add2time("Get k-mer counts from sample");

#Multithreading for checking all subroots
my @newSmallTrees = checkSubroots(@smalltrees);
# Searching from each checked subroot
print STDERR "Searching...\n";
if (scalar @newSmallTrees){
	foreach(@newSmallTrees){
		print STDERR "\n\n---> Subroot $_\n" if $verbose;
		$subroot_constant = subrootConstant($_);
		$current_subroot = $_;
		printStrains($_, $_);
		#print "\n" if $verbose;
	}
}
else {
	print "No match found\n";
}
add2time("Search process");

################################
 ########### CLEANUP ##########
################################
# Remove all unnecessary files #
system "rm $gc_seekerTemp";


add2time("Cleanup");

################################
 ##### COVER CALC/OUTPUT #####
################################
my $covsum = 0; # SUM OF ALL COVERAGES FOR LATER PERCENTAGE
foreach (keys %finalout){
	$covsum+=$finalout{$_}[2];
}

if ($covsum == 0 && scalar keys %finalout != 0){
	push (@warnings, "WARNING: covsum in COVER CALC/OUTPUT  is 0, giving value 1, coverages are not trustworthy");
	$covsum = 1;
}

open(OUTPUT, '>', $outputfilename);
print OUTPUT "Sample:$projectname\n";
foreach my $outBlock (sort { $finalout{$b}[2] <=> $finalout{$a}[2] } keys %finalout) {
   	#@{$finalout{$_}} 0 - known/unknown, 1 - strains hash (strain: found k-mer%), 2 - coverage
	my $percent = sprintf("%.5f",$finalout{$outBlock}[2]/$covsum*100);
	print OUTPUT "$percent%\t$finalout{$outBlock}[0]\t";

	my @strainarray = split(/\n/, $finalout{$outBlock}[1]);
	my %sortinghash;
	foreach (@strainarray) {
		chomp $_;
		my @strain = split(/\t/, $_);
		$sortinghash{$strain[0]} = $strain[1];
	}

	my $i = 1;
	foreach my $str (sort {$sortinghash{$b} <=> $sortinghash{$a} } keys %sortinghash) {
		#my $found = inSample($str);
		if ($i == scalar keys %sortinghash){
			print OUTPUT "$str"; #\tfound $found of $treeinf{$str}[4]\n";

		} else {
			print OUTPUT "$str,"; #\tfound $found of $treeinf{$str}[4]\n";
		}

		$i++;
	}
	print OUTPUT "\n";
}

################################
 ######### PRINT TIME #########
################################
print STDERR "------\n";
my $totalTime = sec2readable(time - $^T);
foreach (@timeArray){
	my ($what, $dura) = split("\t", $_);
	$dura = sec2readable($dura);
	print STDERR "$what: $dura\n";
}
print STDERR "Database root-leaf count: $n_rootleaf\n";
print STDERR "TOTAL RUN TIME: $totalTime\n";
print STDERR "------\n";
print STDERR "Input files:\n".(join("\n", @sample_fasta));
print STDERR "\nResults saved in: $outputfilename\n";

# Print warnings
my $warnNr = scalar @warnings;
if ($warnNr){
	print "----\nWARNING MESSAGES PRODUCED:\n";
	foreach (@warnings){
		print "$_\n";
	}
}
close OUTPUT;
