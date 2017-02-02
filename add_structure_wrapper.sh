#!/bin/bash
# coding: utf8

directory=`dirname $0`

VAR=`echo $3`
VAR0=`echo ${#VAR[0]}`

if [ ${VAR:VAR0-3:VAR0} == 'dat' ]	#Si fichier, modification extension
then
	VAR1=`echo ${VAR:0:VAR0-3}`
	VAR2=`echo $VAR1''out`
	
	mv $3 $VAR2
	
	params=$2" ""--pdb="$4" "$VAR2
	python $directory/add_structure.py $params
	
	fullfilename=$(basename $3)
	
	mv $VAR2 $3	#Remet extension de départ


else
	params=$2" "$3" ""--pdb="$4
	python $directory/add_structure.py $params

	fullfilename=$(basename $4)

fi


#Pour renommer le nom des répertoires et dossiers

extension=${fullfilename##*.}
filename=${fullfilename%.*}

VAR3=`head -1 $4`
VAR4=`echo $VAR3 | sed 's/.* //'`
VAR5=`echo $VAR4`

mv $directory/structures/$2/$filename $directory/structures/$2/$VAR5
mv $directory/structures/$2/$VAR5/$filename.info $directory/structures/$2/$VAR5/$VAR5.info
mv $directory/structures/$2/$VAR5/$filename''_i.dat $directory/structures/$2/$VAR5/$VAR5''_i.dat
mv $directory/structures/$2/$VAR5/$filename''_s.dat $directory/structures/$2/$VAR5/$VAR5''_s.dat

# mettre à jour le nom dans les listes !!
$directory/replacestring.sh $filename $VAR5 $directory/structures/list.dat 
$directory/replacestring.sh $filename $VAR5 $directory/structures/NP_distrib.dat
$directory/replacestring.sh $filename $VAR5 $directory/structures/ABC_s_distrib.dat
$directory/replacestring.sh $filename $VAR5 $directory/structures/ABC_i_distrib.dat
$directory/replacestring.sh $filename $VAR5 $directory/structures/ABC_old_i_distrib.dat
$directory/replacestring.sh $filename $VAR5 $directory/structures/ABC_old_s_distrib.dat




# ---------------------------------------

# Ajout nom protéine dans .loc

sh $directory/makeloc.sh $directory/structures/list.dat ../../../../../tool-data

