#!/usr/bin/python

import sys
import os.path

OUT=open('toto.txt','w')
OUT.write('Protein:' + sys.argv[1] + '\n')
# for the list of structures, we need to take away the name of the protein
# each structures is in the form PROT:STRUC,
# so split each element by : and take the second
strucs=sys.argv[2].split(",")
OUT.write('Structures:')
if len(strucs)>1:
        for s in strucs[:-1]:
                OUT.write(s.split(":")[1]+",")
OUT.write(strucs[-1].split(":")[1])
OUT.write("\n")

print sys.argv

if len(sys.argv) == 3:
	OUT.write('DNA parameter:' + 'ABC_s' + '\n')
	OUT.write('Symmetrization:' + 'No' + '\n')
	OUT.write('Keep tmp:' + 'No' + '\n')
	OUT.write('Temperature factor:' + '1.0' + '\n')

else:
	OUT.write('DNA parameter:' + sys.argv[7] + '\n')
	OUT.write('Symmetrization:' + sys.argv[4] + '\n')
	OUT.write('Keep tmp:' + sys.argv[5] + '\n')
	OUT.write('Temperature factor:' + sys.argv[6] + '\n')
	if sys.argv[3] != '0' :
		OUT.write('Pattern:' + sys.argv[3] + '\n')

OUT.close()

