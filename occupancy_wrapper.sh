#!/bin/bash
# coding: utf8

directory=`dirname $0`
params="--output="$1" "$2" -t "$3

python $directory/occupancy.py $params