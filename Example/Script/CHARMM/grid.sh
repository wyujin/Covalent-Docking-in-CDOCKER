#!/bin/bash/

## parameters for the position and size of the grid
space=0.5
x_center=`cat xcen`
y_center=`cat ycen`
z_center=`cat zcen`
max_len=`cat maxlen`
rcta=`awk '{print $1}' cov_parameter`
rctb=`awk '{print $2}' cov_parameter`
hmax=`awk '{print $3}' cov_parameter`

## parameters for hard core potential
emax=100
mine=-100
maxe=100
$charmm space=$space xcen=$x_center ycen=$y_center zcen=$z_center xmax=$max_len\
 emax=$emax mine=$mine maxe=$maxe rcta=$rcta rctb=$rctb hmax=$hmax < genGrid.inp > genGrid.out

## parameters for hard core potential -- SA
emax=3
mine=-20
maxe=40
$charmm space=$space xcen=$x_center ycen=$y_center zcen=$z_center xmax=$max_len\
 emax=$emax mine=$mine maxe=$maxe rcta=$rcta rctb=$rctb hmax=$hmax < genGrid.inp > genGrid.out

## parameters for hard core potential -- SA
emax=15
mine=-120
maxe=-2
$charmm space=$space xcen=$x_center ycen=$y_center zcen=$z_center xmax=$max_len\
 emax=$emax mine=$mine maxe=$maxe rcta=$rcta rctb=$rctb hmax=$hmax < genGrid.inp > genGrid.out

