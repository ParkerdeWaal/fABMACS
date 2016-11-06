#!/bin/bash
# getReference.sh - Parker de Waal 2016
# extracts reference positions from .gro file for specified atoms
# usage: bash getReference.sh npt.gro
#

# set your CV atoms here in the order in which they will be patched
array=( 2096 2098 2104 2102 2100 2106 2107 2108 )
for i in "${array[@]}"
do
        awk -v id="$i" '{if($3 == id){print $4,$5,$6}}' $1
done
