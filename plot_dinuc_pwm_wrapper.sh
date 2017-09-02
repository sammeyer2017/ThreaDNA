#!/bin/bash
# coding: utf8

directory=`dirname $0`
params="--output="$1" "$2

python $directory/plot_dinuc_pwm.py $params
