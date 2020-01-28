# clade specific prot models
Clade specific models of protein sequence evolution

Files in this repository support:

Pandey, A. and Braun, E.L. (2020) Protein evolution is structure dependent and non-homogeneous 
across the tree of life. Submitted

These models were generated to understand different patterns of protein sequence evolution
in various clades. The models can be used in phylogenetic analyses in IQ-TREE and other
programs.

################################################################################
There are four sets of files in this repository:

I. PAML format amino acid exchangeability matrices

Bird1REV.dat
Bird2REV.dat
MamREV.dat
PlantREV.dat
OomyceteREV.dat
YeastREV.dat
EukREV.dat

These are text files that can be used as empirical models of protein evolution in IQ-TREE or
PAML. With the exception of "EukREV.dat", which was trained on diverse eukaryotes, all of
these models were trained using arbitrarily chosen proteins from specific clades (see "Details 
regarding the training data" below).

--------------------------------------------------------------------------------
II. Nexus file with "models" block containing the XB mixture models

XBmodels_PandeyBraun.nex

The XB (eXposed-Buried) mixture models described by Pandey and Braun (2020) can be implemented
in IQ-TREE using this file. Simply call IQ-TREE with the options:

	-mdef XBmodels_PandeyBraun.nex -m <MODELNAME>
	
	MODELNAMES:

     JarXB   = Bird 1 model
     PrumXB  = Bird 2 model
     MamXB   = Mammal model
     OomyXB  = Oomycete model
     PlantXB = Plant model
     YeastXB = Yeast model
     EukXB   = The "all Euk" model
     EX2llg  = Identical to EX2 from Le et al. (2008), included for comparison
     
Le,S.Q. et al. (2008) Phylogenetic mixture models for proteins. Philos. Trans. Royal Soc. B, 363, 3965-3976.

--------------------------------------------------------------------------------
III. Excel files with amino acid exchangeabilities and other information

exchange_Pandey_Braun.xlsx -- 
	Exchangeabilities for the "all sites" models, including normalized values and comparisons to
	the symmetric EX matrix and matrices for changes in amino acid side-chain volume and polarity.
	The minimum number of nucleotide substitutions necessary for various substitutions is also
	provided in this spreadsheet.

XB_exchange_Pandey_Braun.xlsx --
	Similar to the "exchange_Pandey_Braun.xlsx" filed, but this file contains information about
	the two sub-models of the XB mixture models.

aa_freq_Pandey_Braun.xlsx --
	Equilibrium amino acid frequencies estimated from each training set. This spreadsheet also
	compares the differences between the frequencies for each XB sub-model and the "all sites" 
	frequencies to amino acid side chain polarity. This was done to show that 

EX_matrix_sym.xlsx --
	Symmetric version of the Yampolsky and Stoltzfus (2005) experimental exchangeability matrix.

NOTE: all of these excel spreadsheet include a "README" sheet or a text box with additional information.

Yampolsky,L.Y. and Stoltzfus,A. (2005) The exchangeability of amino acids in proteins. Genetics, 170, 1459-1472.

--------------------------------------------------------------------------------
IV. Programs

fit_XBmixture_models.pl --
	Perl program to identify the best-fitting XB mixture model for individual proteins; the code
	includes extensive comments.

concatenate_nex.py --
	Python program used to concatenate nexus files. Requires biopython
	
NOTE: we used the following pipeline to assign protein structures
	https://github.com/aakanksha12/Structural_class_assignment_pipeline
	
More details about the pipeline can be found in:

Pandey,A. and Braun,E.L. (2019) Phylogenetic analyses of sites in different protein structural environments result in distinct placements of the Metazo-an root. Preprints, 2019100302 http://doi.org/10.20944/preprints201910.0302.v1

################################################################################
Details regarding the training data

All models were optimized using GTR20+FO+I+G4

****
Bird1REV -- Training data:

Input data: 48 sequences with 161122 amino-acid sites
Number of constant sites: 99726 (= 61.8947% of all sites)
Number of invariant (constant or ambiguous constant) sites: 99726 (= 61.8947% of all sites)
Number of parsimony informative sites: 34784
Number of distinct site patterns: 84348

250 proteins selected from:

Jarvis,E.D. et al. (2014) Whole-genome analyses resolve early branches in the tree of life of modern birds. Science, 346, 1320-1331.
Jarvis,E.D. et al. (2015) Phylogenomic analyses data of the avian phylo-genomics project. GigaScience 4:4.

****
Bird2REV -- Training data:

Input data: 317 sequences with 111038 amino-acid sites
Number of constant sites: 61445 (= 55.3369% of all sites)
Number of invariant (constant or ambiguous constant) sites: 61445 (= 55.3369% of all sites)
Number of parsimony informative sites: 38485
Number of distinct site patterns: 70975

All data from Prum et al. (2015), expanded to 317 taxa extracted from genome sequences using the script described in Reddy et al. (2017)

Prum,R.O. et al. (2015) A comprehensive phylogeny of birds (Aves) using targeted next-generation DNA sequencing. Nature, 526, 569-573.
Reddy,S. et al. (2017) Why do phylogenomic data sets yield conflicting trees? Data type influences the avian tree of life more than taxon sampling. Syst. Biol., 66, 857-879.

****
MamREV -- Training data:

Input data: 116 sequences with 238319 amino-acid sites
Number of constant sites: 114423 (= 48.0125% of all sites)
Number of invariant (constant or ambiguous constant) sites: 114423 (= 48.0125% of all sites)
Number of parsimony informative sites: 92889
Number of distinct site patterns: 189543

250 proteins selected from:

Douzery,E.J. et al. (2014) OrthoMaM v8: A database of orthologous exons and coding sequences for comparative genomics in mammals. Mol. Biol. Evol., 31, 1923-1928.

****
PlantREV -- Training data:

Input data: 46 sequences with 80315 amino-acid sites
Number of constant sites: 26374 (= 32.8382% of all sites)
Number of invariant (constant or ambiguous constant) sites: 26374 (= 32.8382% of all sites)
Number of parsimony informative sites: 40822
Number of distinct site patterns: 61058

All data from:

Xi,Z. et al. (2014) Coalescent versus concatenation methods and the placement of Amborella as sister to water lilies. Syst. Biol., 63, 919-932.

****
OomyceteREV -- Training data:

Input data: 17 sequences with 83312 amino-acid sites
Number of constant sites: 38417 (= 46.1122% of all sites)
Number of invariant (constant or ambiguous constant) sites: 38417 (= 46.1122% of all sites)
Number of parsimony informative sites: 33979
Number of distinct site patterns: 35739

Dataset 1, described on page p. 202 of:

Ascunce,M.S. et al. (2017) Phylogenomic analysis supports multiple instances of polyphyly in the oomycete peronosporalean lineage. Mol. Phylogenet. Evol., 114, 199-211.

****
YeastREV -- Training data:

Input data: 343 sequences with 81802 amino-acid sites
Number of constant sites: 9509 (= 11.6244% of all sites)
Number of invariant (constant or ambiguous constant) sites: 9509 (= 11.6244% of all sites)
Number of parsimony informative sites: 68457
Number of distinct site patterns: 76181

277 proteins selected from:

Shen,X.X. et al. (2018) Tempo and mode of genome evolution in the budding yeast subphylum. Cell, 175, 1533-1545.

****
EukREV -- Training data:

Input data: 149 sequences with 58469 amino-acid sites
Number of constant sites: 7415 (= 12.6819% of all sites)
Number of invariant (constant or ambiguous constant) sites: 7415 (= 12.6819% of all sites)
Number of parsimony informative sites: 46893
Number of distinct site patterns: 57004

All data from:

Strassert,J.F. et al. (2019) New phylogenomic analysis of the enigmatic phy-lum Telonemia further resolves the eukaryote tree of life. Mol. Biol. Evol., 36, 757-765.
