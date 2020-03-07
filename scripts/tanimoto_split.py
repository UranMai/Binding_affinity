import argparse
from argparse import ArgumentParser
import sys
sys.path.append('/usr/local/lib/python3.7/site-packages/')

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
import pandas as pd


def data_preparation(path_to_file): #input: csv file
  data = pd.read_csv(path_to_file)
  data = pd.DataFrame(data, columns=['smi', 'Score']) #two interesting columns: smi: 'Smiles' and Score: 'binding affinity'
  data['mol'] = data['smi'].apply(lambda x: Chem.MolFromSmiles(x)) #third column of mol formats
  data = data.dropna()
  data = data.reset_index(drop=True)
  data.loc[data.Score > -30, 'score'] = 0 # the binary scores
  data.loc[data.Score < -30, 'score'] = 1
  return data


def Cluster_Fps(data, cutoff=0.4):
  fps = [AllChem.GetMorganFingerprintAsBitVect(x,2,1024) for x in data['mol']]
  dists = []
  len_fps = len(fps)
  for i in range(1, len_fps):
    sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i]) #Tanimoto similarities(float values) between fps
    dists.extend([1-x for x in sims])

  clusters = Butina.ClusterData(dists, len_fps, cutoff, isDistData = True)
  return clusters


def make_data(cluster, data):
  datum = data.copy()
  data_list = list(sum(cluster, ()))
  length_train = int(len(data_list)*.8)
  length_valid = int(len(data_list)*.9)

  data_train = data_list[:length_train]
  data_valid = data_list[length_train:length_valid]
  data_test = data_list[length_valid:]

  dframe_train = pd.DataFrame()
  dframe_valid = pd.DataFrame() #one parameter in it, column=['smi','Score']
  dframe_test  = pd.DataFrame()

  
  for id, value in enumerate(data_train):
    dframe_train = dframe_train.append(datum.iloc[[value], [0,1,2,3]], ignore_index=True)

  for value in data_valid:
    dframe_valid = dframe_valid.append(datum.iloc[[value], [0,1,2,3]], ignore_index=True)

  for value in data_test:
    dframe_test = dframe_test.append(datum.iloc[[value], [0,1,2,3]], ignore_index=True)
  
  dframe_train.to_csv('train_data.csv')
  dframe_valid.to_csv('valid_data.csv')
  dframe_test.to_csv('test_data.csv')
  
if __name__ == "__main__":
  parser = ArgumentParser()
  parser.add_argument('--data_path', type=str, required=True, help='Path to whole data with smiles and scores')
  parser.add_argument('--cutoff', type=float, default=0.4, help='Cutoff to distinguish clusters')
  args = parser.parse_args()

  data = data_preparation(args.data_path)
  clusters = Cluster_Fps(data, args.cutoff)
  make_data(clusters, data)
