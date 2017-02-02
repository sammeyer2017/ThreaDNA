#!/bin/bash
# Script that replaces a given string in a file by another one, and outputs it either to the same file, or a new one
# Args : oldstring newstring filename (newfilename)

# cat $3 | awk '{ gsub(/'$1'/,/'$2'/);print}' >$4
case "$#" in 
"3" )
cat $3 | awk -v repl=$2 '{ gsub(/'$1'/,repl);print}' >$3.new
mv $3.new $3
;;
"4" )
cat $3 | awk -v repl=$2 '{ gsub(/'$1'/,repl);print}' >$4
;;
esac


exit