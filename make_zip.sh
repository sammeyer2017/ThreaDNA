#!/bin/bash
rm -rf threadna_install
mkdir threadna_install
cp -r add_structure.py calculate_energy.py DNApp.py pwm.py plot_dinuc_pwm.py plot_pwm_with_biopython.py exec.sh globvar.py INSTALL README README.md find_python_version.py seqmotifs.py example occupancy.py struc_list.py util.py param structures threadna_install
chmod u+x threadna_install/INSTALL
chmod u+x threadna_install/exec.sh
chmod u+x threadna_install/find_python_version.py
zip -r threadna_install.zip threadna_install
rm -rf threadna_install
