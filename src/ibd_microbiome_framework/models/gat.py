"""Small, auditable No-Prior GAT for the shared microbiome contract."""
from __future__ import annotations
import torch
from torch import nn
from torch_geometric.nn import GATv2Conv, global_mean_pool

class NoPriorGAT(nn.Module):
    def __init__(self, in_channels: int = 1, hidden: int = 16, heads: int = 2):
        super().__init__()
        self.conv1 = GATv2Conv(in_channels, hidden, heads=heads, dropout=0.1)
        self.conv2 = GATv2Conv(hidden * heads, hidden, heads=1, dropout=0.1)
        self.head = nn.Sequential(nn.Linear(hidden, 16), nn.ReLU(), nn.Dropout(0.1), nn.Linear(16, 1))
    def forward(self, x, edge_index, batch):
        x = torch.relu(self.conv1(x, edge_index))
        x = torch.relu(self.conv2(x, edge_index))
        return self.head(global_mean_pool(x, batch)).view(-1)
