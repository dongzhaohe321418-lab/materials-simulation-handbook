"""CGCNN implementation in PyTorch Geometric.

Extracted from docs/ch10-gnn/03-cgcnn.md (S 10.3.3).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    from torch_geometric.data import Data, Dataset
    try:
        from torch_geometric.loader import DataLoader
    except ImportError:  # older PyG
        from torch_geometric.data import DataLoader  # type: ignore
    from torch_geometric.nn import MessagePassing, global_mean_pool
    _TORCH_AVAILABLE = True
except ImportError:
    _TORCH_AVAILABLE = False


if _TORCH_AVAILABLE:

    class GaussianBasis(nn.Module):
        """Expand a scalar distance r in a fixed Gaussian basis."""

        def __init__(
            self,
            r_min: float = 0.0,
            r_max: float = 8.0,
            n_basis: int = 64,
            sigma: float | None = None,
        ) -> None:
            super().__init__()
            centres = torch.linspace(r_min, r_max, n_basis)
            self.register_buffer("centres", centres)
            if sigma is None:
                sigma = float(centres[1] - centres[0])
            self.sigma = sigma

        def forward(self, r: "torch.Tensor") -> "torch.Tensor":
            delta = r.unsqueeze(-1) - self.centres
            return torch.exp(-0.5 * (delta / self.sigma) ** 2)

    class CGCNNConv(MessagePassing):
        """One CGCNN message-passing layer with edge gating."""

        def __init__(self, atom_dim: int, edge_dim: int) -> None:
            super().__init__(aggr="add")
            z_dim = 2 * atom_dim + edge_dim
            self.gate_linear = nn.Linear(z_dim, atom_dim)
            self.core_linear = nn.Linear(z_dim, atom_dim)
            self.bn_msg = nn.BatchNorm1d(atom_dim)
            self.bn_out = nn.BatchNorm1d(atom_dim)

        def forward(self, h, edge_index, e):
            agg = self.propagate(edge_index, h=h, e=e)
            return self.bn_out(h + agg)

        def message(self, h_i, h_j, e):
            z = torch.cat([h_i, h_j, e], dim=-1)
            gate = torch.sigmoid(self.gate_linear(z))
            core = F.softplus(self.core_linear(z))
            return self.bn_msg(gate * core)

    class CGCNN(nn.Module):
        """Crystal Graph Convolutional Neural Network (Xie & Grossman 2018)."""

        def __init__(
            self,
            n_elements: int = 100,
            atom_dim: int = 64,
            edge_dim: int = 64,
            n_conv: int = 3,
            hidden_dim: int = 128,
            n_targets: int = 1,
        ) -> None:
            super().__init__()
            self.embedding = nn.Embedding(n_elements, atom_dim)
            self.convs = nn.ModuleList(
                [CGCNNConv(atom_dim, edge_dim) for _ in range(n_conv)]
            )
            self.head = nn.Sequential(
                nn.Linear(atom_dim, hidden_dim),
                nn.Softplus(),
                nn.Linear(hidden_dim, n_targets),
            )

        def forward(self, data):
            h = self.embedding(data.Z)
            for conv in self.convs:
                h = conv(h, data.edge_index, data.edge_attr)
            h_G = global_mean_pool(h, data.batch)
            return self.head(h_G).squeeze(-1)


def torch_available() -> bool:
    return _TORCH_AVAILABLE
