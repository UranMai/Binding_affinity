import os
import sys
import argparse
import numpy
import time 

#
#"convert receptor.pdb to pdbqt without HETATM section"
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

parser = argparse.ArgumentParser(description='Protein preparation')
parser.add_argument('-r', '--receptor', nargs='+')
args = parser.parse_args()

def main():
	if args.receptor:
		receptor(sys.argv[2])

if __name__ == '__main__':
	start = time.time()
	main()
	end = time.time()
	print(50*'-')
	print("Execution_time {:.2f} sec".format(end-start))




