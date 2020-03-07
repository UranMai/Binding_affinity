import argparse
from argparse import ArgumentParser
import sys
sys.path.append('/usr/local/lib/python3.7/site-packages/')

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
import pandas as pd


def data_preparation(path_to_file): #input: 'csv file'
  data = pd.read_csv(path_to_file)
  data = pd.DataFrame(data, columns=['smi', 'Score']) #two interesting columns: smi: 'Smiles' and Score: 'binding affinity'
  data['mol'] = data['smi'].apply(lambda x: Chem.MolFromSmiles(x)) #third column of mol formats
  data = data.dropna()
  data = data.reset_index(drop=True)
  data.loc[data.Score > -30, 'score'] = 0 # the binary scores
  data.loc[data.Score < -30, 'score'] = 1
  return data


def generate_scaffold(mol):
    mol = Chem.MolFromSmiles(mol)
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    return scaffold

def train_valid_test_split(data,
          frac_train=.8,
          frac_valid=.1,
          frac_test=.1):
    """
        Splits internal compounds into train/validation/test by scaffold.
        Return the ids  of molecules' scaffold
    """
    
    scaffolds = {}
    for ind, smiles in enumerate(data.smi):
      scaffold = generate_scaffold(smiles)
      if scaffold not in scaffolds:
        scaffolds[scaffold] = [ind]
      else:
        scaffolds[scaffold].append(ind)

    # Sort from largest to smallest scaffold sets (with more same scaffolds,
    scaffolds = {key: sorted(value) for key, value in scaffolds.items()}
    scaffold_sets = [
        scaffold_set
        for (scaffold, scaffold_set) in sorted(
            scaffolds.items(), key=lambda x: (len(x[1]), x[1][0]), reverse=True) #sort by len of sets 
    ]

    train_cutoff = frac_train * len(data)
    valid_cutoff = (frac_train + frac_valid) * len(data)
    train_inds, valid_inds, test_inds = [], [], []
   
    for scaffold_set in scaffold_sets:
      if len(train_inds) + len(scaffold_set) > train_cutoff:
        if len(train_inds) + len(valid_inds) + len(scaffold_set) > valid_cutoff:
          test_inds += scaffold_set
        else:
          valid_inds += scaffold_set
      else:
        train_inds += scaffold_set

    return train_inds, valid_inds, test_inds


def make_data(data):

  train_inds = train_valid_test_split(data)[0]
  valid_inds = train_valid_test_split(data)[1]
  test_inds = train_valid_test_split(data)[2]
  #one parameter in it, column=['smi','Score']
  dframe_train = pd.DataFrame()
  dframe_valid = pd.DataFrame() 
  dframe_test  = pd.DataFrame()

  
  for id, value in enumerate(train_inds):
    dframe_train = dframe_train.append(data.iloc[[value], [0,1,2,3]], ignore_index=True)

  for value in valid_inds:
    dframe_valid = dframe_valid.append(data.iloc[[value], [0,1,2,3]], ignore_index=True)

  for value in test_inds:
    dframe_test = dframe_test.append(data.iloc[[value], [0,1,2,3]], ignore_index=True)

  dframe_train.to_csv('train_data.csv')
  dframe_valid.to_csv('valid_data.csv')
  dframe_test.to_csv('test_data.csv')
  
  

if __name__ == "__main__":
  parser = ArgumentParser()
  parser.add_argument('--data_path', type=str, required=True, help='Path to whole data with smiles and scores')
  args = parser.parse_args()
  data = data_preparation(args.data_path)
  make_data(data)
  
  

