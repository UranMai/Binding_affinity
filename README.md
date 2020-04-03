# Binding_affinity
Drug design, Machine learning, MPNN

# Installing libraries
* Install Anaconda
* Create conda environment
```
  * $ conda create --name <env_name> python={3.6 or 3.7}
  * $ conda activate <env_name>
```
* RDKit ```conda install -q -y -c conda-forge rdkit```
* Pytorch+cpu ```conda install pytorch torchvision cpuonly -c pytorch```
* Pytorch geometric
  ```
  python3 -m pip install torch-scatter==latest+cu101 -f https://pytorch-geometric.com/whl/torch-1.4.0.html
  python3 -m pip install torch-sparse==latest+cu101 -f https://pytorch-geometric.com/whl/torch-1.4.0.html
  python3 -m pip install torch-cluster==latest+cu101 -f https://pytorch-geometric.com/whl/torch-1.4.0.html
  python3 -m pip install torch-spline-conv==latest+cu101 -f https://pytorch-geometric.com/whl/torch-1.4.0.html
  python3 -m pip install torch-geometric
  ```
# Execution
```
python data_prep.py -d data.csv -s sampling -t train_size
python model.py @arguments.txt >> result.txt
```
