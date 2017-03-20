#!/bin/bash

OUTPUTDIR=$PWD/FH_BH_Geometry
INPUTFILE=$PWD/sipm_geometry_tanmay.csv

if [ ! -d $OUTPUTDIR ]; then
    mkdir -p $OUTPUTDIR
fi

# for i in {9..12}; do
#     cat $INPUTFILE | grep "FH${i} " > $OUTPUTDIR/geometry_FH$i.txt
# done

# for i in {1..12}; do
#     cat $INPUTFILE | grep "BH${i} " > $OUTPUTDIR/geometry_BH$i.txt
# done

OFFSET=12
for i in `seq -3 -1`; do
    cat sipm_feb28.txt | grep -e "^${i} " > $OUTPUTDIR/geometry_FH$(($i+$OFFSET)).txt
done

cat sipm_feb28.txt | grep -e "^ 0 " > $OUTPUTDIR/geometry_FH$(($OFFSET)).txt

for i in `seq 1 9`; do
    cat sipm_feb28.txt | grep -e "^ ${i} " > $OUTPUTDIR/geometry_BH${i}.txt
done

for i in `seq 10 12`; do
    cat sipm_feb28.txt | grep -e "^${i} " > $OUTPUTDIR/geometry_BH${i}.txt
done
