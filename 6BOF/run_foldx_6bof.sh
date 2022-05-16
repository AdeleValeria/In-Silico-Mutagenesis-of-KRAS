#!/bin/bash
for f in individual_list*
do
   sub="${f#*-}"
   sub2="${sub%.*}"
   if [[ "$f" == *"hetero"* ]]
   then
      ./foldx5 --command=BuildModel --pdb=6bof_hetero_${sub2}.pdb --mutant-file=$f &&
      ./foldx5 --command=AnalyseComplex --pdb=6bof_hetero_${sub2}_1.pdb --analyseComplexChains=A,B
   else
      ./foldx5 --command=BuildModel --pdb=6bof_${sub2}.pdb --mutant-file=$f &&
      ./foldx5 --command=AnalyseComplex --pdb=6bof_${sub2}_1.pdb --analyseComplexChains=A,B
   fi
done

