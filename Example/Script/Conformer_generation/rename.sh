#!/bin/bash

## This is used to prepare ligand conformer pdb

echo now process conformer $1
sed s/HETATM/ATOM"  "/g tmp.pdb | convpdb.pl -segnames -setchain A | sed s/PRO[0-9A-Z]/LIGA/g | \
	sed s/UNL/LIG/g > tmpligand.pdb 

num=`grep ATOM ligandrtf | sed /LPH/d | wc -l | awk '{print $1}'`
grep ATOM ligandrtf | sed /LPH/d | awk '{print $2}' > tmplist

rm -f tmptmpligand.pdb
for j in `seq 1 $num`; do
	new=`head -n $j tmplist | tail -n 1`
	head -n $j tmpligand.pdb | tail -n 1 | awk '{$3=name; print}' name=$new >> tmptmpligand.pdb
	done

awk '{printf "%-6s %4s %4s %3s %1s %3s %11s %7s %7s %5s %5s %9s \n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' tmptmpligand.pdb > tmpligand.pdb
echo TER >> tmpligand.pdb
echo END >> tmpligand.pdb

convpdb.pl -segnames -setchain A tmpligand.pdb | sed s/PRO[0-9A-Z]/LIGA/g > conformer/${1}.pdb

