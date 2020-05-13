# Smina docking
Smina is fork of Autodock Vina

Resources: 

*  https://sourceforge.net/projects/smina/
* https://github.com/sarisabban/Notes/blob/master/AutoDock.py

# Help
OpenBabel, PyMOL, Smina


# Processing
Create receptor.pdbqt for docking
```
python3 receptor.py -r receptor.pdb
```
If we have txt file with smiles of ligands (smiles.txt), convert to separate **Ligands/*.pdbqt files**
```
python3 preparation.py -d smiles.txt -r receptor.pdb
```

Find the binding box coords and write to **conf.txt**
```
../pymol/pymol bound_box.py -b receptor.pdbqt
```
Create log and pdbqt output and run the best affinities to results.txt
```
bash smina.sh
```
