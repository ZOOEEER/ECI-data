***********************************************************
    Readme file for PredictEnzymeInteraction.exe
***********************************************************

Required files:
--------------------
1.  PredictEnzymeInteraction.exe
2.  PredictEnzymeInteraction.ctf
3.  PredictEnzymeInteraction_main.c
4.  PredictEnzymeInteraction_mcc_component_data.c

System requirements:
---------------------------------
The current version of the software has been compiled to run on Windows XP or Windows 7.  To run the software on a computer that does have Matlab version 7 or higher, the MCRInstaller.exe program should be first run to install the necessary components. 


Data files (examples data files can be found in the data directory) :
--------------
Three files are required to make a prediction:
1) File containing the interaction data for the enzymes with acceptors and donors.   This should be a tab delimited text file, with each donor/acceptor molecules in a row and each enzyme in a column.  The first row of the file should give the names of the enzymes and the first column shuold contain the names of the donors/acceptors.    The interaction data should start in the 3rd row and 3rd column and contain the following numerical values:

    Negative = 0/empty
    Positive  = 1
    Unclear  = 2
    Missing  = 3

2) File containing the protein sequences for enzymes in the interaction file in FASTA format. The order of the enzymes in this file does not matter.

3) File containing the protein sequence for the prediction in FASTA fie. Note that is this file contains multiple protein sequences, a prediction will only be made about the first one. 


