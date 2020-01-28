#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;

############################################################################
# fit_XBmixture_models.pl version 1.0
#
# This program fits a set of mixture models to a protein multiple sequence 
# alignment in relaxed phylip format
#
# As written, the model set is limited to the "XB" (exposed-buried) models
# from Pandey and Braun, submitted. The model set is listed in the
# array @mixmodelname and the hash %mixmodelhash. This program assesses the
# fit of all Pandey and Braun XB models and the the EX2 model of Li et al. 
# 2008 (herein called EX2llg, citation below)
#
# Citations: Le, S. Q., Lartillot, N., & Gascuel, O. (2008) Phylogenetic 
#            mixture models for proteins. Philosophical Transactions of the
#            Royal Society B: Biological Sciences, 363(1512), 3965-3976.
#            doi: 10.1098/rstb.2008.0180
#
#            Pandey, A., Braun, E. L. (submitted) Protein evolution is 
#            structure dependent and non-homogeneous across the tree of 
#            life.
############################################################################

############################################################################
# Dependencies:
#    XBmodels_PandeyBraun: nexus format file with a models block that
#	                       contain the mixture models, must be in the
#                          same directory as this program
#
# Program dependencies:
#    IQ-TREE: absolute dependency, used to calculate the likelihood of 
#             the mixture models. Can also be used to estimate the tree
#             using the best-fitting homogeneous model (e.g., JTT, LG, etc.)
#             within this program before testing the mixture models
#    FastME:  used to generate a minimum evolution tree (using the LG
#             model) if a reference tree is unavailable
#    MPBoot:  used to generate a parsimony tree if a reference tree is
#             unavailable
#
# The reference tree can also be generated using the best-fitting empirical
# model in IQ-TREE (this is slow) or a newick treefile name can be passed on
# the command line (if a reference tree is available)
#
# Enter the paths to the programs for these dependencies below
############################################################################

my($iqexec) = "iqtree"; # path to iqtree executable
my($nthreads) = "-nt 2"; # number of threads for iqtree

my($fastMEexec) = "fastme"; # path to fastme executable
my($mpexec) = "mpboot"; # path to mpboot executable

cleanup(); # run the cleanup routine

############################################################################
# Initialize variables
############################################################################

my($progname) = $0;
my($saveIQ) = 1; # 1 to save the .iqtree files, otherwise files are deleted

my($iter);
my($jter);
my($kter);
my($lter);
my($mter);

my($tempch);
my($tempvar);
my @temparray1;
my @temparray2;

if ( @ARGV < 3 || @ARGV > 4 ) {
	print "Usage:\n  \$ $progname <infile> <outfile> <appendfile> <tree>\n";
	print "  infile     = relaxed phylip format protein alignment\n";
	print "  outfile    = prefix for output files\n";
	print "    FILES: - <prefix>.fit.txt is a tab-delimited list of model fit\n";
	print "             statistics (lnL, AIC, AICc, BIC) and the treelength\n";
	print "           - <prefix>.weights.txt is a file with Akaike weights for\n";
	print "             and models sorted by their fit\n";
	print "           - <prefix>.credibleset.txt is a file listing those models\n";
	print "             in the 95% credible set (based on the AICc)\n";
	print "  appendfile = append short summary of best model to file\n";
	print "  tree       = (optional) newick file with tree topology\n";
	print "    NOTE: - pass FASTME to generate tree using FastME\n";
	print "          - pass PARS to generate tree parsimony (using MPBoot)\n";
	print "          - pass NULL or IQ to generate tree using the best empirical\n";
	print "            model in IQ-TREE (omitting a treefile is the equivalent\n";
	print "            passing the term NULL)\n";
	print "exiting...\n";
	exit;
}

my($infile)=$ARGV[0];
my($outfile)=$ARGV[1];
my($appendfile)=$ARGV[2];
my($treefile)="NULL";

system("cp $infile infile.temp");
if ( @ARGV == 4 ) {
	$treefile=$ARGV[3];
	if ( $treefile eq "IQ" || $treefile eq "IQTREE" ) { $treefile = "NULL"; }
	unless ( $treefile eq "FASTME" || $treefile eq "PARS" ) {
		unless ( $treefile eq "NULL" ) {
			unless ( -e "$treefile" ) { # check whetehr treefile exists
				print "ERROR: cannot find tree file $treefile\nexiting...\n\n";
				exit;
			}
			system("cp $treefile infile.tre");
		}
	}
}

############################################################################
# First, identify the best tree using the best "standard" non-mixture
# model (unless a treefile name was passed)
############################################################################
if ( $treefile eq "NULL" ) {
	print "Generating reference tree topology using IQ-TREE...";
	system("$iqexec -s infile.temp -st AA -pre round1.temp -m TEST $nthreads -redo > round1.temp.screen");
	system("mv round1.temp.treefile infile.tre");
	print "done\n\n";
}
elsif ( $treefile eq "FASTME" ) {
	print "Generating reference tree topology using FastME...";
	if ( -e "infile.tre" ) { unlink("infile.tre"); }
	system("$fastMEexec -i infile.temp -o infile.tre -p L -n -s > round1.temp.screen");
	if ( -z "infile.tre" ) {
		print "\n\nWARNING!!!\nFastME failed to generate tree\n\n";
		open (my $ERRF, ">>$appendfile") or die "Could not open file $appendfile for output\n";
		print $ERRF "$infile\tFASTME_err\t";
		print $ERRF "--\t--\t--\t--\t";
		print $ERRF "--\t--\t--\t--\t--\n";
		close($ERRF) or die "Could not close file $appendfile\n";
		exit;
	}
	print "done\n\n";
}
elsif ( $treefile eq "PARS" ) {
	print "Generating reference tree topology using MPBoot...";
	system("$mpexec -s infile.temp > round1.temp.screen");
	system("mv infile.temp.treefile infile.tre");
	print "done\n\n";
	
}
else { print "Using tree from $treefile\n\n"; }

############################################################################
# Second, iterate through the mixture models. The exposed/buried mixture
# models are assumed to be defined in a nexus file with a models block that
# is called "XBmodels_PandeyBraun.nex"
############################################################################

# First, test whether "XBmodels_PandeyBraun.nex" exists
unless ( -e "XBmodels_PandeyBraun.nex" ) {
	print "ERROR: mixture model definition file XBmodels_PandeyBraun.nex does not exist\nexiting...\n\n";
	exit;
}

# Define the mixture model set
# NOTE: the +I and +I+G models are excluded when mixture models are used
#     - @mixmodel type is 1 for new, 0 for standard (currently available in IQ-TREE)
#     - the standard models (EX2 for the mixture models) must be listed first
my @mixmodelname = ( "EX2llg", "EX2llg+G4", "EX2llg+FO", "EX2llg+FO+G4", 
					 "EukXB", "EukXB+G4", "EukXB+FO", "EukXB+FO+G4",
					 "JarXB", "JarXB+G4", "JarXB+FO", "JarXB+FO+G4",
					 "MamXB", "MamXB+G4", "MamXB+FO", "MamXB+FO+G4",
					 "OomyXB", "OomyXB+G4", "OomyXB+FO", "OomyXB+FO+G4",
					 "PlantXB", "PlantXB+G4", "PlantXB+FO", "PlantXB+FO+G4",
					 "PrumXB", "PrumXB+G4", "PrumXB+FO", "PrumXB+FO+G4",
					 "YeastXB", "YeastXB+G4", "YeastXB+FO", "YeastXB+FO+G4" );
my @mixmodeltype = ( 0, 0, 0, 0,
					 1, 1, 1, 1,
					 1, 1, 1, 1,
					 1, 1, 1, 1,
					 1, 1, 1, 1,
					 1, 1, 1, 1,
					 1, 1, 1, 1,
					 1, 1, 1, 1 );
my($mixmodelnum) = 32;
my($bestmixmodel);
my($bestmixAICc);	my($bestmixK);	my($bestmixlnL);
my($beststdmixmodel);
my($beststdmixAICc);	my($beststdmixK);	my($beststdmixlnL);

# Initialize a hash that will hold delta AICc values (all values initialized as 0)
my %mixmodelhash = (
	"EX2llg" => 0,
	"EX2llg+G4" => 0,
	"EX2llg+FO" => 0,
	"EX2llg+FO+G4" => 0,
	
	"EukXB" => 0,
	"EukXB+G4" => 0,
	"EukXB+FO" => 0,
	"EukXB+FO+G4" => 0,
	
	"JarXB" => 0,
	"JarXB+G4" => 0,
	"JarXB+FO" => 0,
	"JarXB+FO+G4" => 0,
	
	"MamXB" => 0,
	"MamXB+G4" => 0,
	"MamXB+FO" => 0,
	"MamXB+FO+G4" => 0,
	
	"OomyXB" => 0,
	"OomyXB+G4" => 0,
	"OomyXB+FO" => 0,
	"OomyXB+FO+G4" => 0,
	
	"PlantXB" => 0,
	"PlantXB+G4" => 0,
	"PlantXB+FO" => 0,
	"PlantXB+FO+G4" => 0,
	
	"PrumXB" => 0,
	"PrumXB+G4" => 0,
	"PrumXB+FO" => 0,
	"PrumXB+FO+G4" => 0,
	
	"YeastXB" => 0,
	"YeastXB+G4" => 0,
	"YeastXB+FO" => 0,
	"YeastXB+FO+G4" => 0,
);


print "Analyzing fit for $mixmodelnum mixture models:\n";

open (my $OUTF, ">$outfile.fit.txt") or die "Could not open file $outfile.fit.txt for output\n";
print $OUTF "model\tlnL\tK\tAIC\tAICc\tBIC\tTreelen\tIntBL\n";

for ($iter=0; $iter<$mixmodelnum; $iter++) {
	# Optimize model parameters - note that output prefix is fixed as round2.temp (see extract_data subroutine)
	system("$iqexec -s infile.temp -st AA -pre round2.temp -mdef XBmodels_PandeyBraun.nex -m $mixmodelname[$iter] -te infile.tre $nthreads -redo > round2.temp.screen");
	# archive the IQtree files as a single concatenated file
	if ( $saveIQ == 1 ) {
		if ( $iter == 0 ) { system("echo \"##############################\" > $infile.iqtree.txt"); }
		else { 
			system("echo >> $infile.iqtree.txt");
			system("echo \"##############################\" >> $infile.iqtree.txt");
		}
		system("echo \"# Model $iter\" >> $infile.iqtree.txt");
		system("echo \"# $mixmodelname[$iter]\" >> $infile.iqtree.txt");
		system("echo \"##############################\" >> $infile.iqtree.txt");
		system("cat round2.temp.iqtree >> $infile.iqtree.txt"); 
	}
	@temparray1 = extract_data();
#	print "\n\n -- @temparray1\n";
	# Output the model information
	print "  $iter";
	print ". $mixmodelname[$iter]\tlnL = $temparray1[0]\n";
	print $OUTF "$mixmodelname[$iter]";
	for ($kter=0; $kter<7; $kter++) {
		print $OUTF "\t$temparray1[$kter]";
	}
	print $OUTF "\n";
	# Populate the hash with AICc values
	$mixmodelhash{"$mixmodelname[$iter]"} = $temparray1[3];
	# Update best model information
	if ( $iter == 0 ) {
		$bestmixmodel = $mixmodelname[0];
		$bestmixlnL = $temparray1[0];
		$bestmixK = $temparray1[1];
		$bestmixAICc = $temparray1[3];
		# also initialize the best "standard" model (standard models must be listed first in the array)
		$beststdmixmodel = $mixmodelname[0];
		$beststdmixlnL = $temparray1[0];
		$beststdmixK = $temparray1[1];
		$beststdmixAICc = $temparray1[3];
	}
	elsif ( $temparray1[3] < $bestmixAICc ) {
		$bestmixmodel = $mixmodelname[$iter];
		$bestmixlnL = $temparray1[0];
		$bestmixK = $temparray1[1];
		$bestmixAICc = $temparray1[3];
		if ( $mixmodeltype[$iter] == 0 ) {
			$beststdmixmodel = $mixmodelname[$iter];
			$beststdmixlnL = $temparray1[0];
			$beststdmixK = $temparray1[1];
			$beststdmixAICc = $temparray1[3];		
		}
	}
}

close($OUTF) or die "Could not close file $outfile.fit.txt\n";

# Populate the model hash with delta AICc values
for ($iter=0; $iter<$mixmodelnum; $iter++) {
	# Subtract the best AICc from each AICc value in the hash (i.e., generate delta AICc)
	$mixmodelhash{"$mixmodelname[$iter]"} = $mixmodelhash{"$mixmodelname[$iter]"} - $bestmixAICc;
}

# Transfer the delta AICc values to sorted arrays, calculated Akaike weights, establish
# the 95% credible set
my @sortedmixmodels;
my @sortedmixAICc;
my @sortedEXPdelta; # exp( -0.5 * deltaAICc )
my($sumEXPdelta) = 0;
my @sortedmixWT; # Akaike weight and cumulative Akaike weight
my @sortedmixcumWT = 0;
my @credibleset;	$credibleset[0] = 1;

$iter = 0;
foreach my $mixname (sort { $mixmodelhash{$a} <=> $mixmodelhash{$b} } keys %mixmodelhash) {
	$sortedmixmodels[$iter] = $mixname;
	$sortedmixAICc[$iter] = $mixmodelhash{$mixname};
	$tempvar = int($sortedmixAICc[$iter] * 10000 + 0.5); # this rounds to four decimal points
	$sortedmixAICc[$iter] = $tempvar / 10000;
	$sortedEXPdelta[$iter] = exp( -0.5 * $sortedmixAICc[$iter] );
	$sumEXPdelta = $sumEXPdelta + $sortedEXPdelta[$iter];
	$iter++;
}

for ($iter=0; $iter<$mixmodelnum; $iter++) {
	$sortedmixWT[$iter] = $sortedEXPdelta[$iter] / $sumEXPdelta;
	if ( $iter == 0 ) { $sortedmixcumWT[0] = $sortedmixWT[0]; }
	else { 
		if ( $sortedmixcumWT[$iter-1] < 0.95 ) { $credibleset[$iter] = 1; }
		else { $credibleset[$iter] = 0; }
		$sortedmixcumWT[$iter] = $sortedmixcumWT[$iter-1] + $sortedmixWT[$iter];
	}
}

open (my $WTF, ">$outfile.weights.txt") or die "Could not open file $outfile.weights.txt for output\n";
print $WTF "model\tdelta_AICc\tweight\tcum_weight\tcredible\n";

open (my $CRF, ">$outfile.credibleset.txt") or die "Could not open file $outfile.credibleset.txt for output\n";
print $CRF "model\tdelta_AICc\tweight\tcum_weight\n";

for ($iter=0; $iter<$mixmodelnum; $iter++) {
	print $WTF "$sortedmixmodels[$iter]\t$sortedmixAICc[$iter]\t$sortedmixWT[$iter]\t$sortedmixcumWT[$iter]\t$credibleset[$iter]\n";
	if ( $credibleset[$iter] == 1 ) {
		print $CRF "$sortedmixmodels[$iter]\t$sortedmixAICc[$iter]\t$sortedmixWT[$iter]\t$sortedmixcumWT[$iter]\n";	
	}
}

close($WTF) or die "Could not close file $outfile.weights.txt\n";
close($CRF) or die "Could not close file $outfile.credibleset.txt\n";


############################################################################
# Finally, append the results to the "append file"
############################################################################

print "\nBest model = $bestmixmodel\n";
print "($bestmixK free parameters; lnL = $bestmixlnL; AICc = $bestmixAICc)\n\n";
print "Best standard model = $beststdmixmodel\n";
print "($beststdmixK free parameters; lnL = $beststdmixlnL; AICc = $beststdmixAICc)\n\n";

# Append results to the "append file"
# Order of the output in the tab-delimited "append file":
# infile treefile bestmodel bestK best_lnL bestAICc best_Akaike_weight ... continued next line
#     ... best_std best_stdK best_std_lnL best_stdAICc AICc_diff best_std_Akaike_weight std_in_cred_set
#   NOTE: "std" refers to the standard model (in this case the only standar model is the
#         Le et al. 2008 EX2 model, called EX2llg in this program). AICc_diff is the
#         difference between the AICc of the best model and the AICc of the best standard 
#         model
#   NOTE: "std_in_cred_set" is 1 if the best standard model is in the 95% credible set and
#         0 if it is not
my($AICcdiff) = $beststdmixAICc - $bestmixAICc;

open (my $APPF, ">>$appendfile") or die "Could not open file $appendfile for output\n";
print $APPF "$infile\t$treefile\t";
print $APPF "$bestmixmodel\t$bestmixK\t$bestmixlnL\t$bestmixAICc\t";
for ($iter=0; $iter<$mixmodelnum; $iter++) {
	if ( $sortedmixmodels[$iter] eq $bestmixmodel ) { print $APPF "$sortedmixWT[$iter]\t"; }
}
print $APPF "$beststdmixmodel\t$beststdmixK\t$beststdmixlnL\t$beststdmixAICc\t$AICcdiff\t";
for ($iter=0; $iter<$mixmodelnum; $iter++) {
	if ( $sortedmixmodels[$iter] eq $beststdmixmodel ) { print $APPF "$sortedmixWT[$iter]\t$credibleset[$iter]\n"; }
}
close($APPF) or die "Could not close file $appendfile\n";

cleanup();

exit;
############################################################################
# End of main program, subroutines follow
############################################################################



############################################################################
# extract data subroutine
# NOTE: assumes data are in a file named round2.temp.iqtree
############################################################################
sub extract_data {

	# Use grep to extract model fit statistics from round2.temp.iqtree
	system("grep \"^Log-likelihood of the tree:\" round2.temp.iqtree > round2.temp.fitstats");
	system("grep \"^Number of free parameters\" round2.temp.iqtree >> round2.temp.fitstats");
	system("grep \"^Akaike information criterion\" round2.temp.iqtree >> round2.temp.fitstats");
	system("grep \"^Corrected Akaike information criterion\" round2.temp.iqtree >> round2.temp.fitstats");
	system("grep \"^Bayesian information criterion\" round2.temp.iqtree >> round2.temp.fitstats");
	system("grep \"^Total tree length\" round2.temp.iqtree >> round2.temp.fitstats");
	system("grep \"^Sum of internal branch lengths:\" round2.temp.iqtree >> round2.temp.fitstats");
	
	# Read the extracted data
	open (my $MODF, "round2.temp.fitstats") or die "Could not open file round2.temp.fitstats for input\n";
	my @modinfo = <$MODF>; # Read the temporary file with the model information
	close($MODF) or die "Could not close file round2.temp.fitstats\n";
	
	my @modline;
	my @modelfitstats;
	
	@modline = split(/\s+/, $modinfo[0]);
	$modelfitstats[0] = $modline[4]; # Log-likelihood of the tree -- 4
	@modline = split(/\s+/, $modinfo[1]);
	$modelfitstats[1] = $modline[8]; # Number of free parameters -- 8
	@modline = split(/\s+/, $modinfo[2]);
	$modelfitstats[2] = $modline[5]; # Akaike information criterion (AIC) -- 5
	@modline = split(/\s+/, $modinfo[3]);
	$modelfitstats[3] = $modline[6]; # Corrected Akaike information criterion (AICc) -- 6
	@modline = split(/\s+/, $modinfo[4]);
	$modelfitstats[4] = $modline[5]; # Bayesian information criterion (BIC) -- 5

	@modline = split(/\s+/, $modinfo[5]);
	$modelfitstats[5] = $modline[7]; # Total tree length (sum of branch lengths) -- 7
	@modline = split(/\s+/, $modinfo[6]);
	$modelfitstats[6] = $modline[5]; # Sum of internal branch lengths -- 5

	return @modelfitstats;
}

############################################################################
# remove temporary files
# This subroutine removes a number of scratch files used by this program
############################################################################
sub cleanup {
	if ( -e "infile.temp" ) { unlink("infile.temp"); }
	if ( -e "infile.tre" ) { unlink("infile.tre"); }
	if ( -e "infile.temp.log" ) { unlink("iinfile.temp.log"); }
	if ( -e "infile.temp.mpboot" ) { unlink("infile.temp.mpboot"); }
	if ( -e "iinfile.temp.parstree" ) { unlink("iinfile.temp.parstree"); }
	if ( -e "infile.temp_fastme_stat.txt" ) { unlink("infile.temp_fastme_stat.txt"); }
	# NOTE: I assume that if round1.temp.screen or round2.temp.screen exist
	#       there will be other round1 or round2 files. If a temp.screen file
	#       does not exist I assume the directory is clean for all round1 and/or
	#       round2 files
	if ( -e "round1.temp.screen" ) { system("rm round1.temp.*"); }
	if ( -e "round2.temp.screen" ) { system("rm round2.temp.*"); }
}
