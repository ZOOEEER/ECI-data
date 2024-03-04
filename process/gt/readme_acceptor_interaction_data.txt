***********************************************************
    Readme file for Predict_Acceptor_Interaction.exe
***********************************************************

Required files:
--------------------
1.  Predict_Acceptor_Interaction.exe
2.  Predict_Acceptor_Interaction.ctf
3.  predict_acceptor_interaction_main.c
4.  predict_acceptor_interaction_mcc_component_data.c

System requirements:
---------------------------------
The current version of the software has been compiled to run on Windows XP or Windows 7.  To run the software on a computer that does have Matlab version 7 or higher, the MCRInstaller.exe program should be first run to install the necessary components. 


Interaction Data file (example can be found in the data directory) :
----------------------------

The file should contain all the interaction data as well as information about each of the acceptors and be a tab delimited text file. Each row of the file should correspond to an acceptor.

The first 22 columns should be the following information about the acceptors:
ID   - the id number for the acceptor
Name - the name of the acceptor
Family - the family number (1=Flavonoid 2=Coumarin 3=Cytokinin 4=Cinnamic acid 5=Benzoate 6=Jasmonic acid 7=Gibberellins 8=Auxin 9=Abscisic acid 10=Sinapic acid  11=Other)
LogP -   LogP value
Area  - Accessible Area
Volume - volume
COOH - 1 = the molecule has a COOH group, 0 = doesn't have one 
Num OH - the number of OH groups
F3-OH -  0= not present,   greater than 0 = group present
F-5OH	-  0= not present,   greater than 0 = group present
F6-OH	-  0= not present,   greater than 0 = group present
F7-OH	-  0= not present,   greater than 0 = group present
F13-OH	-  0= not present,   greater than 0 = group present
F14-OH	-  0= not present,   greater than 0 = group present
Cm6-OH-  0= not present,   greater than 0 = group present
Cm7-OH  -  0= not present,   greater than 0 = group present	
Ck3-N	-  0= not present,   greater than 0 = group present
Ck7-N	-  0= not present,   greater than 0 = group present
Ck-OH	-  0= not present,   greater than 0 = group present
Cn2-OH-  0= not present,   greater than 0 = group present
Cn3-OH	-  0= not present,   greater than 0 = group present
Cn4-OH-  0= not present,   greater than 0 = group present


Then the enzyme interaction should be entered as follows:

    Negative = 0/empty
    Positive  = 1
    Unclear  = 2
    Missing  = 3





