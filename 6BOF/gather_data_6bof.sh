#!/bin/bash
for f in Summary_6bof_*
do
   if [[ "$f" == *"hetero"* ]]
   then   
      sub="${f#*hetero_}"
      sub2="${sub%_1*}"
      num="${sub2%?}"
      aa="${sub2: -1}"
      interaction=$(awk 'NR == 10 {print $6}' $f)
      total=$(awk 'NR == 10 {print $3}' Average_6bof_hetero_$sub2.fxout)
      delta_g_wt=$(awk 'NR == 11 {print $2}' Raw_6bof_hetero_$sub2.fxout)
      delta_g_mut=$(awk 'NR == 10 {print $2}' Raw_6bof_hetero_$sub2.fxout)
      clash1=$(awk 'NR == 10 {print $4}' $f)
      clash2=$(awk 'NR == 10 {print $5}' $f)
      sasa=$(awk 'NR == 16 {print $3}' SASA_hetero_$sub2.txt)
      sasaA=$(awk 'NR == 19 {print $4}' SASA_hetero_$sub2.txt)
      echo "hetero","$num","$aa","$interaction","$delta_g_wt","$delta_g_mut","$total","$sasa","$sasaA","$clash1","$clash2" >> 6bof_data.csv
   else
      sub="${f#*6bof_}"
      sub2="${sub%_1*}"
      num="${sub2%?}"
      aa="${sub2: -1}"
      interaction=$(awk 'NR == 10 {print $6}' $f)
      total=$(awk 'NR == 10 {print $3}' Average_6bof_$sub2.fxout)
      delta_g_wt=$(awk 'NR == 11 {print $2}' Raw_6bof_$sub2.fxout)
      delta_g_mut=$(awk 'NR == 10 {print $2}' Raw_6bof_$sub2.fxout)
      clash1=$(awk 'NR == 10 {print $4}' $f)
      clash2=$(awk 'NR == 10 {print $5}' $f)
      sasa=$(awk 'NR == 16 {print $3}' SASA_$sub2.txt)
      sasaA=$(awk 'NR == 19 {print $4}' SASA_$sub2.txt)
      echo "homo","$num","$aa","$interaction","$delta_g_wt","$delta_g_mut","$total","$sasa","$sasaA","$clash1","$clash2" >> 6bof_data.csv
   fi
done

