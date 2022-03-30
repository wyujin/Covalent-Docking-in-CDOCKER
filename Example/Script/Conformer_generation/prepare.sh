#!/bin/bash

############################################################################################
## This is used to prepare ligand conformers
############################################################################################


############################################################################################
## Generate conformers
############################################################################################

## Use RDKit to generate conformers
babel -imol2 ligand.mol2 -omol ligand.mol
python property.py > tmpresult
free=`awk '{print $1}' tmpresult`
conj=`awk '{print $2}' tmpresult`
numConfs=`sed 1d conformer.csv | awk '{if ($2 == free && $3 == conj) print $4}' free=$free conj=$conj`
sed -i'.original' s/changeNum/$numConfs/g conformer.py
python conformer.py

## Minimization with CHARMM
rm -rf optimized
mkdir -p optimized
x_center=`cat xcen`
y_center=`cat ycen`
z_center=`cat zcen`
num=`ls conformer/* | wc -l | awk '{print $1}'`
$charmm num=$num xcen=$x_center ycen=$y_center zcen=$z_center < mini.inp > output

