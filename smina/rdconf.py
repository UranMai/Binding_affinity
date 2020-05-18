

from rdkit.Chem import AllChem as Chem

from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit import DataStructs

import sys,string
# import sys
# from rdkit import Chem
# from rdkit.Chem import AllChem as Chem
from optparse import OptionParser

#convert smiles to sdf
def getRMS(mol, c1,c2):
    (rms,trans) = Chem.GetAlignmentTransform(mol,mol,c1,c2)
    return rms

parser = OptionParser(usage="Usage: %prog [options] <input>.smi <output>.sdf")
parser.add_option("--maxconfs", dest="maxconfs",action="store",
                  help="maximum number of conformers to generate per a molecule (default 25)", default="25", type="int", metavar="CNT")
parser.add_option("--sample_multiplier", dest="sample",action="store",
                  help="sample N*maxconfs conformers and choose the maxconformers with lowest energy (default 1)", default="1", type="float", metavar="N")
parser.add_option("--seed", dest="seed",action="store",
                  help="random seed (default 9162006)", default="9162006", type="int", metavar="s")
parser.add_option("--rms_threshold", dest="rms",action="store",
                  help="filter based on rms (default 0)", default="0", type="float", metavar="R")
parser.add_option("--energy_window", dest="energy",action="store",
                  help="filter based on energy difference with lowest energy conformer", default="-1", type="float", metavar="E")
parser.add_option("-v","--verbose", dest="verbose",action="store_true",default=False,
                  help="verbose output")


    
(options, args) = parser.parse_args()

if(len(args) < 2):
    parser.error("Need input and output")
    sys.exit(-1)
    
input = args[0]
output = args[1]
smifile = open(input)
if options.verbose:
    print("Generating a maximum of" + options.maxconfs + "per a mol")
    
sdwriter = Chem.SDWriter(output)
if sdwriter is None:
    print("Could not open ".output)
    sys.exit(-1)
    
for line in smifile:
    toks = line.split()
    smi = toks[0]
    name = "".join(toks[1:])    
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        if options.verbose:
            print(smi)
        try:
            Chem.SanitizeMol(mol)
            mol.SetProp("_Name",name);
            cids = Chem.EmbedMultipleConfs(mol, int(options.sample*options.maxconfs),randomSeed=options.seed)
            if options.verbose:
                print(len(cids) + "conformers found")
            cenergy = []            
            for conf in cids:
                #not passing confID only minimizes the first conformer
                converged = not Chem.UFFOptimizeMolecule(mol,confId=conf)
                if options.verbose:
                    print("Convergence of conformer"+conf+converged)
                cenergy.append(Chem.UFFGetMoleculeForceField(mol,confId=conf).CalcEnergy())
            
            sortedcids = sorted(cids,key = lambda cid: cenergy[cid])
            if(options.rms == 0):
                cnt = 0;
                for conf in sortedcids:
                    if(cnt >= options.maxconfs):
                        break
                    if(options.energy < 0) or cenergy[conf]-cenergy[0] <= options.energy:
                        sdwriter.write(mol,conf)
                        cnt+=1 
            else:
                written = {}
                for conf in sortedcids:
                    if len(written) >= options.maxconfs:
                        break
                    #check rmsd
                    passed = True
                    for seenconf in written.iterkeys():
                        rms = getRMS(mol,seenconf,conf) 
                        if(rms < options.rms) or (options.energy > 0 and cenergy[conf]-cenergy[0] > options.energy):
                            passed = False
                            break
                    if(passed):
                        written[conf] = True
                        sdwriter.write(mol,conf)
        except (KeyboardInterrupt, SystemExit):
            raise                
        except:
            print("Exception occurred: "+sys.exc_info()[0])
    else:
        print("ERROR:"+smi)