help=\
"ThreaDNA: utility program to analyze the mechanical \n
deformation of DNA in high-resolution protein-\n
DNA complexes, and compute mechanical energy profile \n 
of protein binding along a given sequence in a simple \n
and efficient way. \n
\n
usage : threadna [-h] [-a] [-c] [-t] [-g] [-h]\n
\n
-h : prints this message\n
-a : launches the threadna module to add structures\n
-c : launches the energy calculation threadna module\n
-t : prints threadna installation folder\n
-u : updates the list of structures in the database\n 
-g : launches threadna GUI for energy calculation\n
-s : launches the SeqMotif helper program for sequence patterns\n
-o : launches the Occupancy helper program for transforming an energy profile into an occupancy/coverage profile\n
\n
For further help on options -a and -c, type threadna [-a] [-c] -h. \n"
if [ -e $1 ]
then
  python "$DIR/DNApp.py"
  exit 0
fi
if [ $1 = "-a" ]
then
  shift
  python "$DIR/add_structure.py" "$@"
  exit 0
elif [ $1 = "-c" ]
then
  shift
  python "$DIR/calculate_energy.py" "$@"
  exit 0
elif [ $1 = "-t" ]
then
  echo $DIR
  exit 0
elif [ $1 = "-g" ]
then
  python "$DIR/DNApp.py"
  exit 0
elif [ $1 = "-u" ]
then
  python "$DIR/struc_list.py"
  exit 0
elif [ $1 = "-s" ]
then
  shift
  python "$DIR/seqmotifs.py" "$@"
  exit 0
elif [ $1 = "-o" ]
then
  shift
  python "$DIR/occupancy.py" "$@"
  exit 0
else
  echo -e $help
  exit 0
fi
