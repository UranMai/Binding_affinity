#!/bin/bash

start=$(date +%s);

for file in ./Ligands/*;
do
	tmp=${file%.pdbqt};		# without 'pdbqt' 
	name="${tmp##*/}"; 		#stay only ZINC01234
	
	#maybe the not valid output
	smina --config conf.txt --receptor receptor.pdbqt --ligand "$file" --out ${name}_out.pdbqt --log $name.log 

	#write the best results of affinity
	awk '/^[-+]+$/{getline;print FILENAME,$0}' $name.log >> temp; 
	
done;
	
sort temp -nk 3 > Results; 
rm temp; 
mkdir logs; 
#error i mv also receptor.pdbqt, recheck
mv *.log *.pdbqt logs

end=$(date +%s);
diff=$(( $end-$start ))
echo "Execution time $diff seconds"
