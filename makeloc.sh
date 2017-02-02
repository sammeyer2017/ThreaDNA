#!/bin/bash

# This script takes the structures/list.dat list of structures, and creates threadna_protein_names.loc and threadna_structure_names.loc in the appropriate directory for Galaxy: tool_data
# Argument 1: list.dat file (with path)
# Argument 2: directory of the .loc files to be saved (without / at the end)

# protein names
tail -n +2 $1 |cut -f 1|sort -u > $2/threadna_protein_names.loc

# structures names
tail -n +2 $1 |awk -v OFS=":" '{print $1,$3}' |sort > $2/threadna_structure_names.loc
