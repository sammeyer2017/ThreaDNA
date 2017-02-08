#!/bin/bash
rm -rf threadna_install
mkdir threadna_install
cp -r add_structure.py calculate_energy.py DNApp.py exec.sh globvar.py INSTALL README seqmotifs.py example occupancy.py struc_list.py util.py param structures threadna_install
chmod u+x threadna_install/INSTALL
chmod u+x threadna_install/exec.sh
zip -r threadna_install.zip threadna_install
rm -rf threadna_install
