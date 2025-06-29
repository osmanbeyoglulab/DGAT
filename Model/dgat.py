import torch
import torch.nn as nn
from torch_geometric.nn import GATConv



class GATEncoder(nn.Module):
    def __init__(self, in_channels, hidden_dim, num_heads=16, dropout=0.4):
        super(GATEncoder, self).__init__()
        # GAT layers
        self.conv1 = GATConv(in_channels, 2048)
        self.bn1 = nn.LayerNorm(2048)
        self.skip1 = nn.Linear(in_channels, 2048)

        # Multi-head feature attention
        self.num_heads = num_heads
        self.att_heads = nn.ModuleList([
            nn.Sequential(
                nn.Linear(2048, 512),  # Reduce feature size per head
                nn.ReLU(),
                nn.Linear(512, 2048),
                nn.Sigmoid()
            )
            for _ in range(num_heads)
        ])
        self.last_feature_att = None  # To store the latest attention weights

        # Remaining layers
        self.skip2 = nn.Linear(2048, 1024)
        self.conv2 = GATConv(2048, 1024)
        self.bn2 = nn.LayerNorm(1024)

        self.skip3 = nn.Linear(1024, hidden_dim)
        self.conv3 = GATConv(1024, hidden_dim)
        self.bn3 = nn.LayerNorm(hidden_dim)

        self.dropout = nn.Dropout(dropout)
        self.relu = nn.LeakyReLU(0.2)

    def forward(self, x, edge_index, edge_weight=None, return_feature_att: bool = False):
        # ———— First GAT Layer + Skip ————
        identity = self.skip1(x)  
        x1 = self.conv1(x, edge_index)  
        x1 = x1 + identity
        x1 = self.bn1(x1)
        x1 = self.relu(x1)
        x1 = self.dropout(x1)

        # ———— Multi-Head Feature Attention ————
        # Apply each attention head independently, then aggregate
        head_outputs = [head(x1) for head in self.att_heads]
        att_weights = torch.stack(head_outputs, dim=0).mean(dim=0)  # [num_nodes, 2048]
        
        # Store attention weights for later inspection
        self.last_feature_att = att_weights
        
        # Apply attention to the feature channels
        x1 = x1 * att_weights
        
        # ———— Remaining Layers ————
        identity = self.skip2(x1)
        x2 = self.conv2(x1, edge_index)
        x2 = x2 + identity
        x2 = self.bn2(x2)
        x2 = self.relu(x2)
        x2 = self.dropout(x2)

        identity = self.skip3(x2)
        x3 = self.conv3(x2, edge_index)
        x3 = x3 + identity
        x3 = self.bn3(x3)
        x3 = self.relu(x3)
        x3 = self.dropout(x3)

        if return_feature_att:
            return x3, att_weights
        else:
            return x3

class Decoder_mRNA(nn.Module):
    def __init__(self, hidden_dim, output_dim, dropout=0):
        super(Decoder_mRNA, self).__init__()
        self.fc1 = nn.Linear(hidden_dim, 512)
        self.bn1 = nn.LayerNorm(512)
        self.relu1 = nn.LeakyReLU(0.2)
        #self.dropout1 = nn.Dropout(dropout)
        self.skip1 = nn.Linear(hidden_dim, 512)
        self.fc2 = nn.Linear(512, 1024 )
        self.bn2 = nn.LayerNorm(1024 )
        self.relu2 = nn.LeakyReLU(0.2)
        #self.dropout2 = nn.Dropout(dropout)
        self.skip2 = nn.Linear(512, 1024)
        self.fc3 = nn.Linear(1024, output_dim)
        self.bn3 = nn.LayerNorm(output_dim)
        self.relu3 = nn.LeakyReLU(0.2)
        #self.dropout3 = nn.Dropout(dropout)
        self.skip3 = nn.Linear(1024 , output_dim)
        self.linear_activation = nn.Identity()
    def forward(self, z):
        skip_out = self.skip1(z)
        x = self.fc1(z)
        x = self.bn1(x)
        x = self.relu1(x)
        #x = self.dropout1(x)
        x = x + skip_out
        skip_out = self.skip2(x)
        x = self.fc2(x)
        x = self.bn2(x)
        x = self.relu2(x)
        #x = self.dropout2(x)
        x = x + skip_out
        skip_out = self.skip3(x)
        x = self.fc3(x)
        x = self.bn3(x)
        x = self.relu3(x)
        #x = self.dropout3(x)
        x = x + skip_out
        x = self.linear_activation(x)
        return x


class Decoder_Protein(nn.Module):
    def __init__(self, hidden_dim, protein_list, dropout=0.4):

        super(Decoder_Protein, self).__init__()
        self.shared_fc1 = nn.Linear(hidden_dim, 512)
        self.shared_bn1 = nn.LayerNorm(512)
        self.shared_relu1 = nn.LeakyReLU(0.2)
        #self.shared_dropout1 = nn.Dropout(dropout)
        
        self.shared_fc2 = nn.Linear(512, 256)
        self.shared_bn2 = nn.LayerNorm(256)
        self.shared_relu2 = nn.LeakyReLU(0.2)
        #self.shared_dropout2 = nn.Dropout(dropout)
        
        # multi-branches
        self.protein_branches = nn.ModuleDict({
            p: nn.Sequential(
                nn.Linear(256, 64),
                nn.LayerNorm(64),
                nn.LeakyReLU(0.2),
                #nn.Dropout(dropout),
                nn.Linear(64, 1)
            ) for p in protein_list
        })
        self.linear_activation = nn.Identity()

    def forward(self, z):
        # z: [n_cells, hidden_dim]
        x = self.shared_fc1(z)
        x = self.shared_bn1(x)
        x = self.shared_relu1(x)
        #x = self.shared_dropout1(x)
        
        x = self.shared_fc2(x)
        x = self.shared_bn2(x)
        x = self.shared_relu2(x)
        #x = self.shared_dropout2(x)
        
        # 
        outputs = []
        for p in self.protein_branches:
            branch_out = self.protein_branches[p](x)
            outputs.append(branch_out)

        out = torch.cat(outputs, dim=1)
        out = self.linear_activation(out)
        return out

