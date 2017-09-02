#!/bin/bash
# coding: utf8

directory=`dirname $0`
echo $@
python $directory/pwm.py "$@"
