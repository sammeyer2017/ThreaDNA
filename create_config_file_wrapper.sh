#!/bin/bash

directory=`dirname $0`
output_galaxy=$1
prot=$2
struct=$3
patt=$4
sym=$5
tmp=$6
tmp_fact=$7
param=$8


python $directory/create_config_file.py $prot $struct $patt $sym $tmp $tmp_fact $param
cp ./toto.txt $output_galaxy

