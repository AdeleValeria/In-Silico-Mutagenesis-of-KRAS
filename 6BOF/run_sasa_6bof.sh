#!/bin/bash
for f in 6bof_*_1.pdb
do
   sub="${f#*_}"
   sub2="${sub%_1*}"
   freesasa $f > "SASA_${sub2}.txt"
done
