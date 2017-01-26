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

## Dependencies
The program requires Python with the NumPy library. It was successfully tested on Linux and MacOS platforms, with versions of NumPy >= 1.6, and is more efficient on versions >=1.7. 
Computational requirements: the computation time and memory load are proportional to the analyzed sequence length and number of protein-bound nucleotides. For a single protein-DNA complex analyzed, the computation typically takes ~ 200 sec and < 1Gb RAM for a 10 Mb sequence. 

## Description of the input/output files
To get help on command line options, use the command -h for each program/subprogram: e.g. threadna -h, threadna -c -h, ...

# Subprogram : add a new DNA-protein to the database
Usage: 
threadna -a [-p PDBfile] [-r REF] prot_name input
- prot_name: name of the protein in the database. When combining alternate structural models for the same proteins, mind to use the same name (including the case)
- input: input of the DNA basepair/step coordinates in the complex: these can be given as (i) either a NDB ID (http://ndbserver.rutgers.edu) or (ii) a file containing the coarse-grained coordinates of the DNA within the complex, obtained from an analysis program/webserver (.out from 3DNA or the webserver http://w3dna.rutgers.edu, or .lis from Curves+) or (iii) the .pdb atomic coordinates of the protein-DNA complex if the software x3dna is installed on the computer and accessible in the path from the current directory.
- PDBfile: PDB file corresponding to the coordinates given in input. This allows threaDNA to extract additional information on the deformation model, in particular the structure resolution.
- REF: Position of the reference nucleotide in the structure (default = central nucleotide). This is useful when combining/comparing models where the protein binds different numbers of DNA bases: you should then set a reference to align the binding profiles, for instance the index of a basepair that contacts a given protein residue. Typically, for the nucleosome it is the dyad basepair.

RECOMMENDATION FOR THE INPUT: 
- if you are analyzing a protein present in the NDB database, use the option (i) BUT CHECK THAT ALL BASEPAIRS ARE PRESENT IN THE LISTS on the server, as there are sometimes missing basepairs, which will result in errors or absurd deformation energies. This happens particularly in the case of extremely distorted DNA oligomers, where it can be difficult to properly define the base reference frames. 
- if you are analyzing a protein not present or corrupted on the NDB server, download the PDB file and let it run on the Web3DNA webserver to get the .out file. AGAIN, CHECK THAT ALL BASEPAIRS WERE PROPERLY ANALYZED AND ARE PRESENT IN THE LIST. 
- If it does not work or you want to include several structures, try installing Curves+ or contact the authors for help. 

# Main program
Computes an energy profile along a sequence
Usage: 
threadna -c sequence_file conf_file
- sequence_file: FASTA format
- conf_file: description of the calculation parameters: each line is made of a keyword and a value, separated by ":". List of keywords below. 

Output: bedgraph file with 5 columns:
- name of the sequence(s) (better use short names to reduce file size)
- position of the protein-bound reference nucleotide
- position of the protein-bound reference nucleotide
- deformation energy (in a.u. or pseudo k_B T)
additional columns:
- starting nucleotide of the bound protein
- ending nucleotide of the bound protein

PARAMETERS: 

Protein: Name of the protein in the database

Structures: List of alternate protein-DNA structural models that will be combined in the calculation, separated by commas (","). Each structure must have been previously incorporated into the database with the subprogram above, with the same protein name (mind the upper/lowercase). 

Pattern: EMPTY (default) or integer or sequence pattern
Describes which part of the protein-DNA complex analyzed will be included in the calculation, and which sequences of the analyzed sequence will be used. 
- EMPTY: the nucleotide deformations of the whole oligomer are taken into account, and the energy is calculated on the whole input sequence
- integer N (< size of the DNA oligomer in the complex analyzed): only the N basepairs/steps around the reference basepair (see above) are considered in the calculation. 
E.g. for the nucleosome structure 1kx5 that incorporates 147 basepairs, setting a value of 81 results in a deformation energy computed on the 81 central basepairs only, i.e. in the central turn of the DNA. 
- Sequence pattern: restricts the calculation to those positions of the input sequence, which match the pattern (in usual ACGT-RYSW-KMBDHV base code). Pattern size must be shorter than protein size. WARNING: the calculation is MUCH more efficient when carried on the entire input sequence. If you want to restrict to a given sequence pattern, you may prefer to run the computation on the whole sequence first, and then to use the SeqMotifs helper program to find the positions with a given pattern. 

Symmetrization: Yes or No (default)
If the protein is two-fold symmetric, both strands of the bound DNA can be considered as equivalent. Since the protein-DNA structures are generally not perfectly symmetric (especially when the DNA sequence is non-palindromic), if the parameter is set to "Yes", the deformations of both DNA strands are analyzed independently in the program, which makes the resulting deformation model symmetric. 

Keep tmp: Yes or No (default)
If No, only keeps the final energy profile obtained from the combination of all analyzed structural models (in arbitrary units). If Yes, an additional directory is created with all individual energy profiles (in pseudo k_B T unit, where the deformation is assumed to be entirely due to the protein, without any thermal noise). Note that their energy scales can be significantly different. 

Temperature factor: float (default 1.0)
Energy scale used in the Boltzmann factor, for the combination of different structures, in units of (k_B*T). For each sequence, a small value will keep only the most favorable deformation state, and a large value will make all alternate deformation states more equivalent. 

Temperature: integer (default 300)
Temperature (in Kelvin units) for the elastic parameters. If different from 300, a correction of the DNA structure/stiffness parameters is applied, taken from (Meyer et al, Biophys J. 2013). Note that this parameter affects only the DNA elasticity. To compare the binding profile a protein at different temperatures, it is crucial to include an estimation of entropy. 

DNA parameter: "NP" or "ABC_s" or "ABC_i" 
Parameter set for the DNA sequence-dependent basepair/step structure/stiffness. 
- "NP": dataset of basepair step deformations from (Balasubramanian, Xu, Olson, Biophys. J. 2009) based on the analysis of protein-DNA crystallographic structures. Sequence dependence: dinucleotides. 
- "ABC_s": dataset of basepair step deformations from (Pasi et al. NAR 2014) based on 1 microsecond molecular dynamics simulations. Sequence dependence: tetranucleotides. 
- "ABC_i": dataset of internal basepair deformations from (Pasi et al. NAR 2014) based on 1 microsecond molecular dynamics simulations. Sequence dependence: trinucleotides. 
Note: for historical reasons, the previous (50 ns) version of the "ABC" datasets can be used with "ABC_s_old" and "ABC_i_old". 

# Helper programs
SeqMotifs
Simple program to find the positions along a given sequence, where the protein binds a given nucleotide at a given position. 

Usage:
threadna -s [-r REFERENCE] fasta_file nucleotide_type position_on_protein protein_length

For instance, for the Fis protein that binds 15 DNA bases, with the basepair 8 taken as the reference (central basepair), and makes a preferential contact with G at position 0 and C at position 15, you would make two successive runs of SeqMotifs: 
threadna -s -r 8 sequence_file.fasta G 1 15
threadna -s -r 8 sequence_file.fasta C 15 15
The resulting tables have exactly the same structure and indexes as the ThreaDNA profile obtained on the same model (using the pattern "15"). The deformation energies associated to the protein bound with G/C at the desired positions can thus be obtained immediately using MatLab/R/Numpy/Excel, and they can be combined, for instance to introduce a free energy penalty when the correct nucleotide is absent. 

Occupancy
Helper program that computes an occupancy and/or coverage profile. 
Occupancy: probability that the protein occupies a given precise position, defined by the reference basepair of the protein-DNA model (e.g. the dyad basepair of a nucleosome). This profile is simply the Boltzmann inverse of the energy profile. 
Coverage: probability that a given base is covered by a protein, which can then occupy different neighbor positions (depending on the protein size along the DNA). This profile is the occupancy profile convoluted by the protein size. 

Usage: 
threadna -o [-out OUTPUT] [-t EFFTEMP] [-cov COVERAGE] input_file

The input file is an energy profile, typically the output of the main ThreaDNA program. Efftemp is the temperature used in the Boltzmann inversion, in k_B T unit (default 0.3). With the option -cov (coverage), you must also provide the protein size along the DNA, e.g. 147 for the nucleosome. 
