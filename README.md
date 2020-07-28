# clade specific prot models
Clade specific models of protein sequence evolution

Files in this repository support:

Pandey, A. and Braun, E.L. (2020) Protein evolution is structure dependent and non-homogeneous 
across the tree of life. Submitted

These models were generated to understand different patterns of protein sequence evolution
in various clades. The models can be used in phylogenetic analyses in IQ-TREE and other
programs.

##################################################
--------------------------------------------------------------------------------
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
III. Information about the models, including amino acid exchangeabilities

Four excel spreadsheets, each of which either includes a "README" sheet or a text box with additional information:

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
	
Two nexus format treefiles:

	clustered_exchange_Pandey_Braun.tre
	clustered_XBexchange_Pandey_Braun.tre

Both treefiles include comments that will be echoed to the screen if they are executed in PAUP. They can also be visualized in programs such as FigTree.

Reference:

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

Pandey,A. and Braun,E.L. (2019) Phylogenetic analyses of sites in different protein structural environments result in distinct placements of the Metazoan root. Preprints, 2019100302 http://doi.org/10.20944/preprints201910.0302.v1

##################################################
--------------------------------------------------------------------------------
Details regarding the training and validation data

All of the data are available from Zenodo:

Brief summary of the training data:
	  
	Clade		# Proteins/Sites	# Taxa	Model	Citation
	-----		----------------	------	-----	--------
	Birds (1)	250/109,969		48	JTT	[26], [27]
	Birds (2)	250/161,112		317	HIVb	[23]
	Mammals		249/238,319		116	HIVb	[28]
	Plants (1)	310/80,315		46	JTT	[29]
	Oomycetes	277/83,312		17	LG	[30]
	Yeasts		200/81,802		343	LG	[31]
	All Euk (1)	248/58,469		149	LG	[32]
	     
Brief summary of the validation data:

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
     
Le
