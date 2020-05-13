import os
import sys
import argparse
import numpy
import time


def download(filename): 
	'''
	convert smi --> pdbqt
	add name of smi, just in case, for further apps
	delete blank lines for appropriate work in pymol
	'''
	with open(filename, 'r') as file:
		for idx, line in enumerate(file):
			name = idx + 1
			#cmd = "obabel -:'{0}' --gen3d -opdb -O {1}.pdb".format(line, name)#oneofvariants
			
			cmd = "obabel -:'%s' --gen3d -opdbqt -O %s.pdbqt" % (line, name)
			os.system(cmd)
			
			file = str(name) + '.pdbqt'
			
			#cmd1 = "sed -i `s/UNNAMED/{0}/ ` {1}.pdb".format(line, name)
			#cmd1 = "sed -i `s/UNNAMED/%s/` %s.pdb" % (line, name)
			#cmd1 = "sed -i 's/UNNAMED/"+line+"/' "+ file
			#print(cmd1)
			#os.system(cmd1)
			
			#name the smiles in pdb files
			with open(file, mode='r') as rf:
				change = rf.read().replace('Name = ', 'Name = '+str(line))
				with open(file, mode='w') as wf:
					wf.write(change)

			#delete blank lines
			cmd1 = "sed -i '/^$/d' %s.pdbqt" % (name) 
			os.system(cmd1)
	os.system("mkdir Ligands")
	os.system("mv *.pdbqt Ligands")
	os.system("cd Ligands, rm receptor.pdbqt")
			
def download1(filename):
	'''
	convert smi --> pdb --> pdbqt
	add name of smi, just in case, for further apps
	delete blank lines for appropriate work in pymol 
	'''
	with open(filename, 'r') as file:
		for idx, line in enumerate(file):
			name = idx + 1
			cmd = "obabel -:'%s' --gen3d -opdb -O %s.pdb" % (line, name)
			###os.system("obabel -:'$line' --gen3d -opdb -O tes.pdb")
			os.system(cmd)
			
			file = str(name) + '.pdb'
						
			#name the smiles in pdb files
			with open(file, mode='r') as rf:
				change = rf.read().replace('UNNAMED', line)
				with open(file, mode='w') as wf:
					wf.write(change)		
			
			#delete blank lines
			cmd1 = "sed -i '/^$/d' %s.pdb" % (name) 
			os.system(cmd1)
			
			cmd2 = "obabel %s -O %s.pdbqt" % (file, name)
			os.system(cmd2)
	os.system("mkdir Ligands")
	os.system("mv *.pdbqt Ligands")

#
#convert receptor.pdb to pdbqt without HETATM
#

def receptor(file):
	os.system("grep ATOM %s > R.pdb" % file)
	cmd = "echo 'load R.pdb; remove resn HOH; h_add elem O or elem N; save protein.pdb; quit' > pymol.pml"
	os.system(cmd)
	os.system("../pymol/pymol pymol.pml")
	os.system("obabel protein.pdb -xr -O temp.pdbqt")
	os.system("grep ATOM temp.pdbqt > receptor.pdbqt")
	cmd1 = "rm R.pdb; rm temp.pdbqt; rm protein.pdb; rm pymol.pml"
	os.system(cmd1)

#parser.add_argument('-d', '--download', type=str, nargs='+', help='Work with smiles')
#parser.add_argument('-r', '--receptor', type=str, nargs='+', help='Process receptor')
#nargs doesnot work with argparse, only with sys.argv

parser = argparse.ArgumentParser(description='Preparation full data')
parser.add_argument('--download', '-d', type=str, help='Prepare smiles ligands')
parser.add_argument('--receptor', '-r', type=str, help='Prepare receptor')
args = parser.parse_args()

def main():
	download(args.download)
	receptor(args.receptor)


if __name__ == '__main__':
	start = time.time()
	main()
	end = time.time()
	print(50*'-')
	print("Execution_time {:.2f} sec".format(end-start))