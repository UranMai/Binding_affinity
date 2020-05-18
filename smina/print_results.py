import os
import sys
import argparse
import numpy
import time 
import pandas as pd

def main(file):
  with open(file, 'r') as f:
    lines = f.readlines()
    scores = []
    first_compound = lines[0].rstrip('\n')
    comps = [first_compound]
    for idx, line in enumerate(lines):
      if '> <minimizedAffinity>' in line:
        score = lines[idx+1].rstrip('\n')
        scores.append(score)

    for idx, line in enumerate(lines):
      if '$$$$' in line:
        names = lines[idx+1].rstrip('\n')
        comps.append(names)
    
    data = {'name' : comps, 'scores' : scores}
    df = pd.DataFrame(data)
    result = df.groupby(['name'])['scores'].max()
    output = result.sort_values(ascending=False)
    return output


parser = argparse.ArgumentParser(description='Read file with the results')
parser.add_argument('-f', '--file', type=str)
args = parser.parse_args()

if __name__ == "__main__":
  file1 = open('results.txt', 'w')
  file1.write(str(main(args.file)))
  #print(main(args.file))