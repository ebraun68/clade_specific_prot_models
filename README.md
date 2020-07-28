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

Pandey, A.; Braun, E.L. Phylogenetic Analyses of Sites in Different Protein Structural Environments Result in Distinct Placements of the Metazoan Root. Biology 2020, 9, 64. https://doi.org/10.3390/biology9040064
(also available in preprint for as: Preprints, 2019100302 http://doi.org/10.20944/preprints201910.0302.v1)

##################################################
--------------------------------------------------------------------------------
Details regarding the training and validation data (all of which are available from Zenodo)

Zenodo link: Pandey, Akanksha, & Braun, Edward L. (2020). Protein evolution is structure dependent and non-homogeneous across the tree of life. http://doi.org/10.5281/zenodo.3964471

Brief summary of the training data:
	  
	Clade		# Proteins/Sites	# Taxa	Model	Citation
	-----		----------------	------	-----	--------
	Birds (1)	250/109,969		48	JTT	[1], [2]
	Birds (2)	250/161,112		317	HIVb	[3]
	Mammals		249/238,319		116	HIVb	[4]
	Plants (1)	310/80,315		46	JTT	[5]
	Oomycetes	277/83,312		17	LG	[6]
	Yeasts		200/81,802		343	LG	[7]
	All Euk (1)	248/58,469		149	LG	[8]
	     
Brief summary of the validation data:
	  
	Clade		# Proteins	# Taxa	Citation
	-----		----------	------	--------
	Birds (1)	200		48	[1], [2]
	Mammals		200		116	[4]
	Plants (2)	200		107	[9]
	Oomycetes	150		15	[6]
	Yeasts		200		343	[7]
	All Euk (2)	149		104	[10]
     
Citations for datasets:

1. E. D. Jarvis, S. Mirarab, A. J. Aberer, B. Li, P. Houde, C. Li, S. Y. W. Ho, B. C. Faircloth, B. Nabholz, J. T. Howard, A. Suh, C. C. Weber, R. R. da Fonseca, J. Li, F. Zhang, H. Li, L. Zhou, N. Narula, L. Liu, G. Ganapathy, B. Boussau, Md. S. Bayzid, V. Zavidovych, S. Subramanian, T. Gabaldón, S. Capella-Gutiérrez, J. Huerta-Cepas, B. Rekepalli, K. Munch, M. Schierup, B. Lindow, W. C;.Warren, D. Ray, R. E Green, M. Bruford, X. Zhan, A. Dixon, S. Li, N. Li, Y. Huang, E. P. Derryberry, M. F. Bertelsen, F. Sheldon, R. T. Brumfield, C. Mello, P. V. Lovell, M. Wirthlin, J. A. Samaniego, A. M. V. Velazquez, A. Alfaro-Núñez, P. F. Campos, T. Sicheritz-Ponten, A. Pas, T .Bailey, P. Scofield, M. Bunce, D. Lambert, Q. Zhou, P. Perelman, A. C. Driskell, G. Ruby, B. Shapiro, Z. Xiong, Y. Zeng, S. Liu, Z. Li, B. Liu, K. Wu, J. Xiao, X. Yinqi, Q. Zheng, Y. Zhang, H. Yang, J. Wang, L. Smeds, F. E. Rheindt, M. Braun, J. Fjeldså, L. Orlando, K. Barker, K. A. Jønsson, W. Johnson, K.-P. Koepfli, S. O'Brien, D. Haussler, O. A. Ryder, C. Rahbek, E. Willerslev, G. R. Graves, T. C. Glenn, J. McCormack, D. Burt, H. Ellegren, P. Alström, S. V. Edwards, A. Stamatakis, D. P. Mindell, J. Cracraft, E. L. Braun, T. Warnow, Wang J., M. T. P. Gilbert, and G. Zhang, “Whole-genome analyses resolve early branches in the tree of life of modern birds.,” Science, vol. 346, no. 6215, pp. 1320–1331, Dec. 2014.

2. E. D. Jarvis, S. Mirarab, A. J. Aberer, B. Li, P. Houde, C. Li, S. Y. W. Ho, B. C. Faircloth, B. Nabholz, J. T. Howard, A. Suh, C. C. Weber, R. R. da Fonseca, A. Alfaro-Núñez, N. Narula, L. Liu, D. Burt, H. Ellegren, S. V. Edwards, A. Stamatakis, D. P, Mindell, J. Cracraft, E. L. Braun, T. Warnow, Wang J., M. T. P. Gilbert, G. Zhang, and Avian Phylogenomics Consortium, “Phylogenomic analyses data of the avian phylogenomics project.,” Gigascience, vol. 4, p. 4, Feb. 2015.

3. R. O. Prum, J. S. Berv, A. Dornburg, D. J. Field, J. P. Townsend, E. M. Lemmon, and A. R. Lemmon, “A comprehensive phylogeny of birds (Aves) using targeted next-generation DNA sequencing.,” Nature, vol. 526, no. 7574, pp. 569–573, Oct. 2015.

4. E. J. P. Douzery, C. Scornavacca, J. Romiguier, K. Belkhir, N. Galtier, F. Delsuc, and V. Ranwez, “OrthoMaM v8: a database of orthologous exons and coding sequences for comparative genomics in mammals.,” Mol. Biol. Evol., vol. 31, no. 7, pp. 1923–1928, Jul. 2014.

5. Z. Xi, L. Liu, J. S. Rest, and C. C. Davis, “Coalescent versus concatenation methods and the placement of Amborella as sister to water lilies.,” Syst. Biol., vol. 63, no. 6, pp. 919–932, Nov. 2014.

6. M. S. Ascunce, J. C. Huguet-Tapia, A. Ortiz-Urquiza, N. O. Keyhani, E. L. Braun, and E. M. Goss, “Phylogenomic analysis supports multiple instances of polyphyly in the oomycete peronosporalean lineage.,” Mol. Phylogenet. Evol., vol. 114, pp. 199–211, Jun. 2017.

7. X.-X. Shen, D. A. Opulente, J. Kominek, X. Zhou, J. L. Steenwyk, K. V. Buh, M. A. B. Haase, J. H. Wisecaver, M. Wang, D. T. Doering, J. T. Boudouris, R. M. Schneider, Q. K. Langdon, M. Ohkuma, R. Endoh, M. Takashima, R.-I. Manabe, N. Čadež, D. Libkind, C. A. Rosa, and A. Rokas, “Tempo and mode of genome evolution in the budding yeast subphylum.,” Cell, vol. 175, no. 6, pp. 1533-1545.e20, Nov. 2018.

8. J. F. H. Strassert, M. Jamy, A. P. Mylnikov, D. V. Tikhonenkov, and F. Burki, “New phylogenomic analysis of the enigmatic phylum Telonemia further resolves the eukaryote tree of life.,” Mol. Biol. Evol., vol. 36, no. 4, pp. 757–765, Apr. 2019.

9. N. J. Wickett, S. Mirarab, N. Nguyen, T. Warnow, E. Carpenter, N. Matasci, S. Ayyampalayam, M. S. Barker, J. G. Burleigh, M. A. Gitzendanner, B. R. Ruhfel, E. Wafula, J. P. Der, S. W. Graham, S. Mathews, M. Melkonian, D. E. Soltis, P. S. Soltis, N. W. Miles, C. J. Rothfels, L. Pokorny, A. Jonathan Shaw, L. DeGironimo, D. W. Stevenson, B. Surek, J. C. Villarreal, B. Roure, H. Philippe, C. W. dePamphilis, T. Chen, M. K. Deyholos, R. S. Baucom, T. M. Kutchan, M. M. Augustin, J. Wang, Y. Zhang, Z. Tian, Z. Yan, X. Wu, X. Sun, G. K.-S. Wong, and J. Leebens-Mack, “Phylotranscriptomic analysis of the origin and early diversification of land plants.,” Proc. Natl. Acad. Sci. USA, vol. 111, no. 45, pp. E4859-68, Nov. 2014.

10. G. Lax, Y. Eglit, L. Eme, E. M. Bertrand, A. J. Roger, and A. G. B. Simpson, “Hemimastigophora is a novel supra-kingdom-level lineage of eukaryotes.,” Nature, vol. 564, no. 7736, pp. 410–414, Nov. 2018.
