# For building k-mer tree databases from Newick tree files and bacteria fasta files for StrainSeeker's seeker.pl
# First build: 21jan2015 by M2rt Roosaare and Mihkel Vaher
# IMPORTANT!!! FASTA file names must be exactly the same as the names in the tree file, with suffix .fna.

use strict;
use warnings;
use Getopt::Long;

my $prevTime = time;

# Get StrainSeeker path
my $path;
if (index($0, "/") == -1){
	$path = "./";
} else {
	my @pathArray = split("/",$0);
	pop(@pathArray);
	$path = (join ("/", @pathArray))."/";
}

############################################
 ################# INPUT ##################
############################################

# Files, manditory
my $newick;
my $dir_strains; # Strains directory
my $dir_tree; # Output lists dir, user defined

# Program parameters, optional
my $word = 32; # Word size for tree database
my $min_kmers = 250; # Minimal amount of k-mers required for potential root nodes
my $child_greater_parent = 250; # Number of times child has to be with greater k-mer count to discard current node for subroot
my $threads = 32; # Number of parallel processes
my $max_in_list = 100000; # Maximum numbers of k-mers in list
my $blacklist; # Blacklist is subtracted only when given
my $help;
my $version;
GetOptions(
    'n=s' => \$newick,
    'newick=s' => \$newick,
    'd=s' => \$dir_strains,
    'dir=s' => \$dir_strains,
    'o=s' => \$dir_tree,
    'output=s' => \$dir_tree,
    'w=i' => \$word,
    'word=i' => \$word,
    'm=i' => \$min_kmers,
    'min=i' => \$min_kmers,
    'g=i' => \$child_greater_parent,
    'greater=i' => \$child_greater_parent,
    't=i' => \$threads,
    'threads=i' => \$threads,
    'h' => \$help,
    'help' => \$help,
    'v' => \$version,
    'version' => \$version,
    'max=i' => \$max_in_list,
    'b=s' => \$blacklist,
    'blacklist=s' => \$blacklist,
     ) or die printHelp()."\n";

if ($version){
	die "builder.pl v1.5\n";
}

if ((!$newick) || (!$dir_strains) || (!$dir_tree) || ($help)){
	die printHelp()."\n";
}

############################################
 ############## GLOBAL PARAMS #############
############################################

# Glist-program locations that are used
my $glm = $path."glistmaker"; # GListmaker
my $glc = $path."glistcompare"; # Glistcompare
my $glq = $path."glistquery"; # Glistquery
my $gc  = $path."gmer_counter"; #GmerCounter


# Dir/file locations
my $unionfile = "union_$word\_union.list"; # makeUnion outputfile
my $tree_data_final = "info.txt"; # Final info file for searching, relationships of nodes

# Globals used for nwk parsing
my $idx = 0;
my @parse_leaves;
my @parse_nodes;
my @parse_paths;
my @table; # Table for nodes and leaves, used in subs parse and converTable

#Other globals
my %treehash; # Key: $node; value: array [Current node] [Parent] [Child 1] [Child 2] [COUNT]
my @timeArray;
my %remain; # List of all nodes that are in subroots and will NOT be removed
my @warnings;
my $blacklistUsed = "False";

############################################
 ########### CHECK IF ALL IS OK ###########
############################################

# Check word size
if ($word > 32 || $word < 2){
	die "Word size has to be in range 2-32\n";
}

if ($threads <= 0){
	$threads = 1;
}

# Check if $dir_tree ends with "/"
while (substr($dir_strains, -1) eq "/"){
	chop $dir_strains;
}

# Check if files exist
if (!-e $glm) { die "GListmaker not found!\n" }
if (!-e $glc) { die "GListcompare not found!\n" }
if (!-e $glq) { die "GListquery not found!\n" }
if (!-e $gc) { die "gmer_counter not found!\n" }
if (($blacklist) && (!-e $blacklist)) { die "Blacklist not found!\n" }
if (!-e $newick) { die "Newick file $newick not found!\n" }
if (!-d $dir_strains) { die "Directory $dir_tree with fasta files not found!\n" }

# Check if blacklist is with same word size as planned db
if (($blacklist) && ((split(/\s+/, (qx/$glq $blacklist -stat/)[1]))[1] != $word)){
	die "Blacklist and database have to be with the same word size (k-mer length)!\n"
}

# Try to make dir for tree list files
if (!-e $dir_tree) {
	die "Cannot make directory: mkdir $dir_tree" if system "mkdir $dir_tree";
	die "Cannot give permissions to directory: chmod 750 $dir_tree" if system "chmod 750 $dir_tree";
}
else {die "Project dir $dir_tree already exists!\n"}

# FASTAS ARE CHECKED IN fastasExist sub

###########################################
 ################# SUBS ##################
###########################################
sub add2time { # Argument: timestamp's name to be added
	my $out = "$_[0]\t".(time - $prevTime);
	push (@timeArray, $out);
	$prevTime = time;
}

sub blacklist { # Argument: array of all nodes; compares given node to blacklist, removes found k-mers from node
	print "COMPARING WITH BLACKLIST...\n";
	my @nodes = @_;
	my $cmd;
	my $rm = "rm ";
	my $count = 0;
	my $total = 0;
	my $n_ofall = scalar @nodes;
	foreach my $node (@nodes) {
		$count++;
		$cmd .= "$glc $dir_tree/$node\_$word\.list $blacklist -d -o $dir_tree/$node & ";
		$rm .= "$dir_tree/$node\_$word\.list ";
		#Multithread execution if thread limit reached
		if ($count == $threads) {
			$total = $total + $threads;
			print "Blacklist: $total of $n_ofall\n";
			if(system "$cmd wait"){
				die "Problem with $cmd in blacklist";
			}

			if(system "$rm"){
				die "Problem with $rm in blacklist";
			}
			$count = 0;
			$cmd = "";
			$rm = "rm ";
		}
	}
	# Check for unexecuted commands (last nodes if number of threads is not a multiple of node count)
	if ($cmd ne "") {
		if (system "$cmd wait"){
			die "Problem with cmd: $cmd in blacklist sub";
		}
		if (($rm ne "rm ") && (system "$rm")){
			push (@warnings, "WARNING: problem with removing in blacklist sub");
		}
	}

	print "Renaming all diff files\n";
	foreach my $list (@nodes) {
		if(system "mv $dir_tree/$list\_$word\_0_diff1.list $dir_tree/$list\_$word\.list"){
			die "Problem with renaming mv $dir_tree/$list\_$word\_0_diff1.list $dir_tree/$list\_$word\.list in blacklist";
		}
	}
}

sub convertTable { # Converts $newick, fills %treehash: @([Current node] [Parent] [Child 1] [Child 2] [COUNT(NA)]), returns @leaves, uses parse()
	my @parse_leaves;
	my @parse_nodes;
	my @parse_paths;

	open (my $ifs, $newick) or die $!;
	while (my $line = <$ifs>) {
		chomp($line);
		$line =~ s/\r//;
		parse ($line, 1, "", 0);
	}
	close ($ifs);
	my $out;
	my @out2;
	my @output; #output array
	my %nodes; #hash of all nodes

	foreach (@table) {
		my $check_root = 1;
		my @current = split(/,/, $_); #order: [node][child1][child2]
		$nodes{$current[0]} = 0; #Put all nodes with children inside hash
		foreach (@table) { #Finds parent for current node
			my @tmp = split(/,/, $_);
			if (($tmp[1] eq $current[0]) || ($tmp[2] eq $current[0]) ) {
				$out = "$current[0]\t$tmp[0]\t$current[1]\t$current[2]\n"; #ALL NODE LINES
				@out2 = ($current[0],$tmp[0],$current[1],$current[2],"NA");
				@{$treehash{$current[0]}} = @out2;
				$check_root = 0;
				last;
			}
		}
		if ($check_root == 1) {
			$out = "$current[0]\tNA\t$current[1]\t$current[2]\n";
			@out2 = ($current[0],"NA",$current[1],$current[2],"NA");
			@{$treehash{$current[0]}} = @out2;
		}
	}
	my %leaves; #hash of all endpoint nodes (strains)
	foreach (@table) { #Finds nodes without children (strains) and outputs them
		my @current = split(/,/, $_);
		if (!exists($leaves{$current[1]})) {
			if (!exists($nodes{$current[1]})) { #Not a node and not yet checked leaf - print out
				$out = "$current[1]";
				@out2 = ($current[1],$current[0],"NA","NA","NA");
				@{$treehash{$current[1]}} = @out2;

				push(@output, $out);
				$leaves{$current[1]} = 0;
			}
		}
		if (!exists($leaves{$current[2]})) {
			if (!exists($nodes{$current[2]})) { #Not a node and not yet checked leaf - print out
					$out = "$current[2]";
	 					@out2 = ($current[2],$current[0],"NA","NA","NA");
				@{$treehash{$current[2]}} = @out2;

				push (@output, $out);
				$leaves{$current[2]} = 0;
			}
		}
	}

	return @output;
}

sub eliminate { # Arguments: node and its child; eliminate redundant k-mers from children nodes; used in node list generation
	my $node = shift;
	my $child = shift;
	if(system "$glc $dir_tree/$child\_$word\.list $dir_tree/$node\_$word\.list -d -o $dir_tree/$child"){
		die "Problem with $glc $dir_tree/$child\_$word\.list $dir_tree/$node\_$word\.list -d -o $dir_tree/$child";
	}
	die "Cannot remove old child list file: rm $dir_tree/$child\_$word\.list" if system "rm $dir_tree/$child\_$word\.list"; #Remove old child list file
	die "Cannot rename file: mv $dir_tree/$child\_$word\_0_diff1.list $dir_tree/$child\_$word\.list" if system "mv $dir_tree/$child\_$word\_0_diff1.list $dir_tree/$child\_$word\.list";
}

sub eliminateMulti { # Arguments: union file where are all to compare with, array that has all nodes from a subtree
	print "ELIMINATING K-MERS OCCURING IN MULTIPLE NODES\n";
	my $union = shift;
	my @nodes = @_;
	my $cmd;
	my $rm = "rm ";
	my $count = 0;
	my $total = 0;
	my $n_ofall = scalar @nodes;
	foreach my $node (@nodes) {
		$count++;
		$cmd .= "$glc $dir_tree/$node\_$word\.list $union -du -o $dir_tree/$node & ";
		$rm .= "$dir_tree/$node\_$word\.list ";
		#Multithread execution if thread limit reached
		if ($count == $threads) {
			$total = $total + $threads;
			print "Eliminate multiple occurring: $total of $n_ofall\n";
			if(system "$cmd wait"){
				die "Problem with $cmd in eliminateMulti";
			}

			if(system "$rm"){
				die "Problem with $rm in eliminateMulti";
			}

			$count = 0;
			$cmd = "";
			$rm = "rm ";
		}
	}
	#check for unexecuted commands (last nodes if number of threads is not a multiple of node count)
	if ($cmd ne "") {
		if (system "$cmd wait"){
			die "Problem with cmd: $cmd in eliminateMulti";
		}
		print "Eliminate multiple occurring: all remaining\n";

		if (($rm ne "rm ") && (system "$rm")){
			die "Problem with removing in eliminateMulti";
		}
	}

	foreach (@nodes) {
		die "Cannot rename file: mv $dir_tree/$_\_$word\_0_diff1.list $dir_tree/$_\_$word\.list" if system "mv $dir_tree/$_\_$word\_0_diff1.list $dir_tree/$_\_$word\.list";
	}
}

sub fastasExist { # Arguments: $dir_strains, @leaves
	my $dir_strains = shift;
	my @leaves = @_;
	foreach (@leaves){
		if (-e "$dir_strains/$_\_$word\.list"){ # If exists list instead of fasta, prefer this
			system "cp $dir_strains/$_\_$word\.list $dir_tree/$_\_$word\.list";
		}
		elsif (!-e "$dir_strains/$_\.fna"){
			die "$dir_strains/$_\.fna does not exist! Wrong directory or .fna file names do not match with .nwk names\n";
		}
	}
}

sub getChildren { # Argument: node; return array of nodes including given itself
	(my $current) = @_;
	my $child1 = $treehash{$current}[2];
	my $child2 = $treehash{$current}[3];
	my @children;
	push (@children, $current);
	if ($child1 eq "NA"){ # if current is strain
		return @children;
	} else {
		push (@children, (getChildren($child1),getChildren($child2)));
		return @children;
	}
}

sub getKmerCount { # Argument: node, get unique k-mer count from a list
	my $node = shift;
	$node = "$dir_tree/".$node."\_$word\.list";
	my @kmer_count = qx/$glq $node -stat/;

	foreach (@kmer_count) {
		if ($_ =~ /Unique/) {
			chomp $_;
			my @tmp = split(/\s+/, $_);
			return $tmp[1];
		}
	}
}

sub genNodes { # Argument: @leaves; make union of children to create parent and remove parent's k-mers from children
	print "Creating node lists: union of children and remove redundant...\n";
	my @nodes = @_; # At first contains only leaves
	my %checked;
	my %forCmd;
	my $cmd;
	my $count = 0;
	my $totalnodes = (scalar keys %treehash);
	my $tot = 0;
	while (@nodes){
		$tot++;
		my $node = shift @nodes;
		my $parent = $treehash{$node}[1];

		if ((!-e "$dir_tree/$node\_$word\.list") && ($treehash{$node}[2] eq "NA")){
			die "$dir_tree/$node\_$word\.list list not found! (leaves not made?)";
		} elsif (($parent ne "NA") && (-e "$dir_tree/$treehash{$parent}[3]\_$word\.list") && (-e "$dir_tree/$treehash{$parent}[2]\_$word\.list") && ($treehash{$node}[2] eq "NA")){
			if (!exists($checked{$parent})) {
				push (@nodes,$parent);
				$checked{$parent} = 1;
			}
			next;
		} elsif ($treehash{$node}[2] ne "NA") {
			if ($node eq "NA"){die "Somehow NA got into building queue, aborting"}
			$forCmd{$node} = 1;
			$count++;
		}

		if ((!@nodes) || ($count == $threads)){
			foreach my $toList (keys %forCmd){
				$cmd .= "$glc $dir_tree/$treehash{$toList}[2]\_$word\.list $dir_tree/$treehash{$toList}[3]\_$word\.list -i -o $dir_tree/$toList & ";
			}

			print "Generating nodes: $tot of $totalnodes\n";
			if(system "$cmd wait"){
				die "Problem with genNodes main command: $cmd";
			}

			foreach (keys %forCmd){
				if(system "mv $dir_tree/$_\_$word\_intrsec.list $dir_tree/$_\_$word\.list"){
					die "Problem with renaming mv $dir_tree/$_\_$word\_intrsec.list $dir_tree/$_\_$word\.list\nlast command $cmd";
				}
				eliminate("$_", "$treehash{$_}[2]"); #Child 1 	#Eliminate redundant k-mers from children nodes and add child kmer precount (not final) to hash
				eliminate("$_", "$treehash{$_}[3]"); #Child 2
				if (($treehash{$_}[1] ne "NA") && (-e "$dir_tree/$treehash{$treehash{$_}[1]}[3]\_$word\.list") && (-e "$dir_tree/$treehash{$treehash{$_}[1]}[2]\_$word\.list") && (!exists($checked{$treehash{$_}[1]}))){
					push (@nodes,$treehash{$_}[1]);
					$checked{$treehash{$_}[1]} = 1;
				}
			}
			undef %forCmd;
			$count = 0;
			$cmd = "";
		}
	}
}

sub leaves2lists { # Argument: @leaves, make lists from all leaves
	print "CREATING LISTS FOR LEAVES...\n";
	my @leaves = @_;
	my $cmd;
	my $count = 0;
	my $total = 0;
	my $len = scalar @leaves;

	foreach my $node (@leaves) {
		if(-e "$dir_tree/$node\_$word\.list"){
			$total++;
			next;
		}

		$cmd .= "$glm $dir_strains/$node\.fna -w $word -o $dir_tree/$node & ";
		$count++;
		#Multithread execution if thread limit reached
		if ($count == $threads) {
			$total = $total + $threads;
			print "Generating leaves: $total of $len\n";
			if(system "$cmd wait"){
				die "Leaves2lists: problem with $cmd";
			}

			$count = 0;
			$cmd = "";
		}
	}
	# Check for unexecuted commands (last nodes if number of threads is not a multiple of node count)
	if ($cmd && $cmd ne "") {
		die "Problem with cmd: $cmd wait" if system "$cmd wait";
	}
}

sub makeUnion { # Argument: array of nodes to make union of
	my @nodes = @_;
	my %forCmd;
	my $cmd;
	my $rm = "rm ";
	my $left = (scalar @nodes)-1;
	my $unioncounter = 0;
	my $last_union;

	if (scalar @nodes == 1){
		die "Problem with cmd: cp $dir_tree/$nodes[0]\_$word\.list union\_$word\_union.list" if system "cp $dir_tree/$nodes[0]\_$word\.list union\_$word\_union.list";
		return;
	} elsif(scalar @nodes == 0){
		return;
	}

	while (@nodes){
		if (scalar @nodes > 1){
			my $unionname = "MULTITHREADED_union_".$unioncounter;
			my $node1 = shift @nodes;
			if (!exists($treehash{$node1})){
				$node1.= "_$word\_union.list";
			} else {
				$node1.= "_$word.list";
			}

			my $node2 = shift @nodes;
			if (!exists($treehash{$node2})){
				$node2.= "_$word\_union.list";
			} else {
				$node2.= "_$word.list";
			}
			@{$forCmd{$unionname}} = ($node1, $node2);
			$unioncounter++;
		}
		if ((scalar @nodes <= 1) || (scalar keys %forCmd >= $threads)){
			print "Making union of lists, remaining: $left\n";
			foreach (keys %forCmd){
				$cmd .= "$glc $dir_tree/$forCmd{$_}[0] $dir_tree/$forCmd{$_}[1] -u -o $dir_tree/$_ & ";
				if ("$dir_tree/$forCmd{$_}[0]" =~ /MULTITHREADED_union_/) {
					$rm .= "$dir_tree/$forCmd{$_}[0] ";
				}
				if ("$dir_tree/$forCmd{$_}[1]" =~ /MULTITHREADED_union_/) {
					$rm .= "$dir_tree/$forCmd{$_}[1] ";
				}
				push (@nodes, $_);
				$last_union = $_;
			}
			if(system "$cmd wait"){
				die "Problem with main command: $cmd";
			}

			if(($rm ne "rm ") && (system "$rm")){
				die "Problem with main command: $rm";
			}

			$left-= scalar keys %forCmd;
			undef %forCmd;
			$cmd = "";
			$rm = "rm ";
			if (scalar @nodes == 1){
				last;
			}
		}
	}
	if (system "mv $dir_tree/$last_union\_$word\_union.list union\_$word\_union.list"){
				die "Problem with command: mv $dir_tree/$last_union\_$word\_union.list union\_$word\_union.list";
	}
	if (glob("$dir_tree/MULTITHREADED_union_*") && (system "rm $dir_tree/MULTITHREADED_union_*")){
		push (@warnings, "Problem with removing: rm $dir_tree/MULTITHREADED_union_*");
	}
}

sub parse { # Parse node | s is the next index from "(" | returns the position of ")", used by converTable()
	my $str = $_[0];
	my $s = $_[1];
	my $path = $_[2];
	my $level = $_[3];
	my $e;
	my $leaf;
	my $pp;

	# Create id
	my $id = "N$idx";
	$idx += 1;
	push (@parse_nodes, $id);

	# Parse left
	my $leftid;
	if (substr ($str, $s, 1) eq "(") {
		# Left is list
		$leftid = "N$idx";
		$e = parse ($str, $s + 1, "$path$id ",	$level + 1) + 1;
	} else {
		# Left is leaf
		$e = index ($str, ":", $s);
		if ($e >= 0) {
			$leaf = substr ($str, $s, $e - $s);
			push (@parse_leaves, $leaf);
			$pp = "$path$id $leaf";
			push (@parse_paths, $pp);
			$leftid = $leaf;
		} else {
			die "Not valid nwk file\n";
		}
	}
	$s = $e;
	# s is now at colon
	$s += 1;

	# Find comma
	$e = index ($str, ",", $s);
	if ($e >= 0) {
		$s = $e + 1;
	} else {
		die "Not valid nwk file\n";
	}

	# Parse right
	my $rightid;
	if (substr ($str, $s, 1) eq "(") {
		# Right is list
		$rightid = "N$idx";
		$e = parse ($str, $s + 1, "$path$id ", $level + 1) + 1;
	} else {
		# Right is leaf
		$e = index ($str, ":", $s);
		if ($e >= 0) {
			$leaf = substr ($str, $s, $e - $s);
			push (@parse_leaves, $leaf);
			$pp = "$path$id $leaf";
			push (@parse_paths, $pp);
			$rightid = $leaf;
		} else {
			die "Not valid nwk file\n";
		}
	}
	$s = $e;
	# s is now at colon
	$s += 1;

	# Find closing brace
	$e = index ($str, ")", $s);
	if ($e >= 0) {
		$s = $e + 1;
	} else {
		die "Not valid nwk file\n";
	}

	push(@table, "$id,$leftid,$rightid");
	return $s;
}

sub printHelp {
	print "Usage: $0 -n <NWK FILE> -d <DIR OF FASTA FILES> -o <USER DEFINED DB NAME> [OPTIONAL PARAMETERS]\n";
	print "Options:
	-h, --help\t - Print this help
	-v, --version\t - Print version of the program
	-n, --newick\t - Guide tree in newick format (same names as fasta files without suffix .fna)
	-d, --dir\t - Directory of fasta files (.fna)
	-o, --output\t - User defined database name
	-b, --blacklist\t - .list file of k-mers unwanted in database (human, plasmids etc)
	-w, --word\t - K-mer length used in database building and later searching (default $word)
	-m, --min\t - Minimal amout of k-mers in node to be considered as subroot (default $min_kmers)
	-g, --greater\t - Maximum times child could have more k-mers than parent (default $child_greater_parent)
	-t, --threads\t - Number of cores used
	-max\t - Maximum number of k-mers in one list (default $max_in_list)\n";
	return "";
}

sub reduceLists { # Argument: array of nodes, randomly removes k-mers until list has $max_in_list k-mers, counts k-mers and adds to %treehash
	print "REDUCING LISTS...\n";
	my @nodes = @_;
	my $cmd;
	my $count = 0;
	my $total = 0;
	my $n_ofall = scalar @nodes;
	my $alreadyDone = 0;
	my @reduced;
	foreach my $node (@nodes) {
		my $initial_count = getKmerCount($node);
		die "reduceLists: could not get inital count for $node" if (!defined $initial_count);
		if ($initial_count <= $max_in_list){ #if less than $max_in_list - ignore
			$alreadyDone++;
			next;
		}

		$cmd .= "$glc $dir_tree/$node\_$word\.list -ss rand_unique $max_in_list -o $dir_tree/$node & ";
		push (@reduced, $node);
		$count++;

		#Multithread execution if thread limit reached
		if ($count == $threads) {
			$total = $total + $threads + $alreadyDone;
			print "Reduce lists: $total of $n_ofall\n";
			if(system "$cmd wait"){
				die "Problem with $cmd in reduceLists";
			}

			$alreadyDone = 0;
			$count = 0;
			$cmd = "";
		}
	}
	# Check for unexecuted commands (last nodes if number of threads is not a multiple of node count)
	if ($cmd ne "") {
		if (system "$cmd wait"){
			die "Problem with last cmd: $cmd in reduceLists";
		}
	}

	print "Renaming all subset files and get new counts\n";
	foreach (@reduced) {
		# Rename
		if(system "mv $dir_tree/".$_."\_subset_$word\.list $dir_tree/".$_."\_$word\.list"){
			die "Problem in reduceLists command: mv $dir_tree/".$_."\_subset_$word\.list $dir_tree/".$_."\_$word\.list";
		}

		# Get count
		$treehash{$_}[4] = $max_in_list;
	}
}

sub sec2readable { # Argument: time in seconds, returns in mins or hours
		my $secs = shift;
		if	  ($secs >= 365*24*60*60) { return sprintf '%.1fy', $secs/(365*24*60*60) }
		elsif ($secs >=		24*60*60) { return sprintf '%.1fd', $secs/(24*60*60) }
		elsif ($secs >=		   60*60) { return sprintf '%.1fh', $secs/(60*60) }
		elsif ($secs >=			  60) { return sprintf '%.1fm', $secs/(60) }
		else	{ return sprintf '%.1fs', $secs}
}

sub subroot { # Argument: node; returns array of subroots; for finding all the roots with minimal k-mer count cutoff
	my ($current) = @_;
	my $child1 = $treehash{$current}[2];
	my $child2 = $treehash{$current}[3];
	my @roots;
	if ($treehash{$current}[4] > $min_kmers) {
		if ($child1 eq "NA") {
			push (@roots, $current);
			return @roots;
		} elsif (($treehash{$child1}[4] / $child_greater_parent > $treehash{$current}[4]) || ($treehash{$child2}[4] / $child_greater_parent > $treehash{$current}[4])){
			my @c1array = subroot($child1);
			my @c2array = subroot($child2);
			if ($c1array[0]){
				push (@roots, @c1array);
			}
			if ($c2array[0]){
				push (@roots, @c2array);
			}
			return @roots;
		} else {
			push (@roots, $current);
			return @roots;
		}
	} elsif ($child1 eq "NA") {
		print "Subroot: WARNING: $current contains less than $min_kmers (".$treehash{$current}[4].") and will be left out, consider lowering 'min' parameter\n";
		delete $treehash{$current};
		push (@warnings, "Subroot: WARNING: $current contains less than $min_kmers (".$treehash{$current}[4].") and will be left out, consider lowering 'min' parameter");
		return;
	} else {
		my @c1array = subroot($child1);
		my @c2array = subroot($child2);
		if ($c1array[0]){
			push (@roots, @c1array);
		}
		if ($c2array[0]){
			push (@roots, @c2array);
		}
		return @roots;
	}
}


###########################################################
 ##################### Main program ######################
###########################################################

#Convert Newick tree into our tree format
print "CONVERTING NEWICK...\n";
my @leaves = convertTable(); # Fills %treehash, returns @leaves

fastasExist($dir_strains,@leaves);
add2time("Start and nwk convert");

leaves2lists(@leaves);
add2time("Leaves generation");

if($blacklist){ # Compare with blacklist only when it is given
	$blacklistUsed = "True";
	blacklist(@leaves);
	add2time("Comparing with blacklist");
}

genNodes(@leaves); # Works through the tree, getting k-mer counts and eliminating redundant k-mers, NODE LIST GENERATION
if (!-e "$dir_tree/N0_$word\.list"){die "$dir_tree/N0_$word\.list was not made"} # Check if root is made
add2time("Node generation and counts");

# Multiple occurance removal
print "CREATING LARGE UNION FOR MULTIPLE OCCURANCE REMOVAL...\n";

makeUnion(keys %treehash);
add2time("Large union generation for multiple occurance removal");
eliminateMulti($unionfile, keys %treehash);
system "rm $unionfile";

add2time("Multiple occurance removal");

# Get count for all leaves and nodes
print "GET INITAL K-MER COUNT IN LISTS...\n";
foreach (keys %treehash) {
	my $tmpCount = getKmerCount($_);
	die "Count loop: getKmerCount for $_ returned nothing\n" if (!defined $tmpCount);
	$treehash{$_}[4] = $tmpCount;
}

# ///////// NODES GENERATED, KMERS COUNTED, MULTIPLE OCCURRING REMOVED

# Find subroots (if main root is OK, it will be the only subroot)
my @subroots = subroot("N0");
if (!scalar @subroots){
	die "No subroots made (database not made), min parameter ($min_kmers) too high?\n";
}
add2time("Getting subroots");

# Remove unused lists
my $rootleaf = 0;
print "REMOVING UNUSED LISTS...\n";
foreach my $root (@subroots){ # Build hash of what to keep
	$treehash{$root}[1] = "NA"; # Change subroot's parent to NA
	if ($treehash{$root}[2] eq "NA"){
		$rootleaf++;
	}
	foreach my $node (getChildren($root)){
		$remain{$node} = 1;
	}
}
foreach (keys %treehash){
	if (!exists($remain{$_})){
		delete $treehash{$_};
		if (system "rm $dir_tree/$_\_$word\.list"){
			push (@warnings, "Cannot remove redundant file: rm $dir_tree/$_\_$word\.list");
		}
	}
}
add2time("Redundant removal");

# Trim very large lists randomly to given value ($max_in_list) and get kmer count
die "No lists remaining, redundant removal removed too much?\n" if (scalar keys %treehash == 0);
reduceLists(keys %treehash);
add2time("List reducing");

# CREATE INFO FILE
my $totalNodes = scalar keys %treehash;
my $totalSubroots = scalar @subroots;

print "Writing tree data to $tree_data_final\n";
open(TREEINF, '>', "$dir_tree/$tree_data_final");
# MAKE HEADER FOR INFO FILE
print TREEINF "Database:$dir_tree\tWordSize:$word\tSubroot_Kmer_cutoff:$min_kmers\tNodes:$totalNodes\tSubroots:$totalSubroots\tRoot-leaves:$rootleaf\tBlacklist_used:$blacklistUsed\tMax_kmers_in_any_list:$max_in_list\tMin_kmers_in_subroots:$min_kmers\tTimes_child_greater_than_parent_to_discard_as_subroot:$child_greater_parent\tNewick_file_used:$newick\tInput_fasta_dir:$dir_strains\n---\n";
# FILL INFO FILE: <current> <parent> <child1> <child2> <count>
foreach my $key (keys %treehash){
	my $out = join("\t",@{$treehash{$key}});
	print TREEINF "$out\n";
}

print TREEINF "---"; # Separate tree relations and subroots

foreach (@subroots) { # Append subroots to file
	print TREEINF "\n$_";
}

# CONVERT DATABASE TO BINARY
print "CREATING BINARY FILE...\n";
# Convert to suitable input format
my $suffix = localtime($^T);
$suffix =~ s/[\s:]+//g;
my $tempFile = "StrainSeeker_temporary_text_file".$suffix.".txt";

print "Writing to temp text file $tempFile\n";
opendir(D, $dir_tree) || die "Can't open directory $dir_tree: $!\n";
open(my $fh, '>>', $tempFile) or die "Could not open file '$tempFile' $!";
while (my $f = readdir(D)) {
    next if (index($f, ".") == -1);

    my ($name, $suffix) = (split(/\./, $f));

    # Remove k-mer len suffix (_32)
    my @jupid = split(/\_/, $name);
    pop(@jupid);
    $name = join("_",@jupid);
    $f = $dir_tree."/".$f;

    next if ($suffix ne "list");

    my $count =  ((split (/\n/, (`$glq $f -stat`)))[2]);
    $count = (split (/\t/, $count))[1];
    print "$f\n";
    print $fh "$name\t$count\t";
    system ("$glq $f | cut -f1 | head -n -2 | tr '\n' '\t'>> $tempFile");

    print $fh "\n";

}
closedir(D);
close $fh;

# Remove redundant
system ("rm $dir_tree/*.list");

# Convert to binary
print "Converting into binary file...\n";
system ("$gc -db $tempFile -w $dir_tree/db_binary");
system ("rm $tempFile");
add2time("Binary db creation");

# TIME STUFF
print "--------------------\n";
my $totalTime = sec2readable(time - $^T);
foreach (@timeArray){
	my ($what, $dura) = split("\t", $_);
	$dura = sec2readable($dura);
	print "$what: $dura\n";
}
print "Total run time: $totalTime\n";

# Print warnings
my $warnNr = scalar @warnings;
if ($warnNr){
	print "----\nWARNING MESSAGES PRODUCED:\n";
	foreach (@warnings){
		print "$_\n";
	}
}

# End
print "--------------------\nDONE\nTotal $totalNodes nodes of which $totalSubroots are subroots of which $rootleaf are root-leaves.\nWarning messages produced: $warnNr\n";
close TREEINF;