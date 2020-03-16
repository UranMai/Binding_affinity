"""import torch
print(torch.__version__)
torch.cuda.is_available()
from torch.nn import Sequential as Seq, Linear, ReLU
from torch_geometric.nn import MessagePassing
from torch_geometric.utils import remove_self_loops
from torch_geometric.utils import add_self_loops
from torch_geometric.data import Data, DataLoader
import mpnn.py
from mpnn.py import mol2graph

import torch
import torch.nn.functional as F
from torch.nn import Sequential, Linear, ReLU, GRU

import torch_geometric.transforms as T
from torch_geometric.nn import NNConv, Set2Set #NNConv (implemented MPNN); Set2Set 
from torch_geometric.data import DataLoader
from torch_geometric.utils import remove_self_loops"""




class Data_prep(object):
	train = pd.read_csv('train_data.csv', index_col=0)
	valid = pd.read_csv('valid_data.csv', index_col=0)
	test = pd.read_csv('test_data.csv', index_col=0)

	train_data = train.drop(train[train.Score>0].index)
	valid_data = valid.drop(valid[valid.Score>0].index)
	test_data = test.drop(test[test.Score>0].index)

	train_data['mol'] = train_data['smi'].apply(lambda x: Chem.MolFromSmiles(x))
	train_data = train_data.dropna()
	train_data = train_data.sample(frac=1) # for shuffling
	
	valid_data['mol'] = valid_data['smi'].apply(lambda x: Chem.MolFromSmiles(x))
	valid_data = valid_data.dropna()
	valid_data = valid_data.sample(frac=1)

	test_data['mol'] = test_data['smi'].apply(lambda x: Chem.MolFromSmiles(x))
	test_data = test_data.dropna()
	test_data = test_data.sample(frac=1)

	train_X = [mol2graph.mol2vec(Chem.MolFromSmiles(m)) for m in train_data.smi]
	for i, data in enumerate(train_X):
	    y = train_data.Score.values[i]
	    data.y = torch.tensor([y], dtype=torch.float)

	valid_X = [mol2graph.mol2vec(Chem.MolFromSmiles(m)) for m in valid_data.smi]
	for i, data in enumerate(valid_X):
	    y = valid_data.Score.values[i]
	    data.y = torch.tensor([y], dtype=torch.float)

	test_X = [mol2graph.mol2vec(Chem.MolFromSmiles(m)) for m in test_data.smi]
	for i, data in enumerate(test_X):
	    y = test_data.Score.values[i]
	    data.y = torch.tensor([y], dtype=torch.float)

	train_loader = DataLoader(train_X, batch_size=128, shuffle=True)
	valid_loader = DataLoader(valid_X, batch_size=128, shuffle=False)
	test_loader = DataLoader(test_X, batch_size=128, shuffle=False)	

class Net(torch.nn.Module):
  def __init__(self, node_dim=75, dim=64):
    super(Net, self).__init__()
    
    self.lin0 = torch.nn.Linear(node_dim, dim) 
    nn = Sequential(Linear(6, 128),  
                    ReLU(), 
                    Linear(128, dim*dim)) 
    self.conv = NNConv(in_channels=dim, out_channels=dim, nn=nn, aggr='add')
    self.gru = GRU(dim, dim)

    self.set2set = Set2Set(dim, processing_steps=3, num_layers=1) 
    self.lin1 = torch.nn.Linear(2*dim, dim)
    self.drop = torch.nn.Dropout(p=0.2)
    self.lin2 = torch.nn.Linear(dim, 1)

  def forward(self, data): 
    out = F.relu(self.lin0(data.x)) 
    out = self.drop(out)
    h = out.unsqueeze(0) 

    for i in range(5):
            m = F.relu(self.conv(out, data.edge_index, data.edge_attr)) 
            out, h = self.gru(m.unsqueeze(0), h)
            out = out.squeeze(0) 

    out = self.set2set(out, data.batch) 
    out = F.relu(self.lin1(out)) 
    out = self.drop(out)
    out = self.lin2(out)

    return out.view(-1) 


device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = Net().to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min',
                                                       factor=0.8, patience=5,
                                                       min_lr=0.00001)

def train(epoch):
    model.train()
    loss_all = 0

    for data in train_loader:
        data = data.to(device)
        optimizer.zero_grad()
        loss = F.mse_loss(model(data), data.y)
        loss.backward()
        loss_all += loss.item() * data.num_graphs
        optimizer.step()
    return loss_all / len(train_loader.dataset)


def test(loader):
    model.eval()
    error = 0
    
    for data in loader:
        data = data.to(device)
        pred = model(data)
        error += F.l1_loss(pred, data.y).item()
    return error / len(loader.dataset)

def target(loader):
	model.eval()
	targets = dict()
	for data in loader:
	  data = data.to(device)
	  y_pred = model(data)
	  for i in range(len(data.y)):
	      targets[data.y[i].item()] = y_pred[i].tolist()
	return targets

if __name__ == '__main__':
	main()	
