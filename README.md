# ThreaDNA  version 1.0

## Installation
If the correct dependencies are installed (see below), you only need to execute the INSTALL script (either from the command line or by double-clicking). 
An executable file "threadna" should appear in your folder, which you can execute by double-clicking. 
If any of these files should open as a text file (with script code) instead of being executed, try to righ-click on the file, and change the file properties to make it executable by the user. 
In the graphical user interface, help is displayed for all options. 

## Principle of the program
The program calculates the DNA deformation energy associated to the binding of a protein, along an input sequence. In the typical case, the protein binds 10-20 basepairs, and the energy profile is computed along a genomic sequence. 

The main program takes the sequence of interest (in the FASTA format) and a configuration file describing the protein-DNA complex model and a few computation parameters. 
IMPORTANT: the model is built from one or several high-resolution structures of the DNA-protein complex of interest: before any calculation can be carried out, these structures must first be incorporated into the program database with the dedicated subprogram (see below). This step must done only once for each structure. 

The calculation is based on the hypothesis that the protein imposes DNA's shape at the nanometer-scale. The bending energy is then calculated from the deformed configuration, using the sequence-dependent structure/elastic parameters for the internal base-pair or the base-pair step deformations. For proteins where several alternate structures are available (typically obtained with different DNA sequences), these can be combined in the calculation: for a given sequence, DNA then "chooses" the weight of each conformation from the different Boltzmann weights. 
