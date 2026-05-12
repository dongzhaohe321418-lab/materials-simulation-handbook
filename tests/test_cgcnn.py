"""Test the CGCNN model trains and reduces loss on a tiny synthetic batch."""
from __future__ import annotations

import pytest

torch = pytest.importorskip("torch")
pyg = pytest.importorskip("torch_geometric")

from tier1.ch10.cgcnn import torch_available


@pytest.mark.skipif(not torch_available(), reason="PyTorch / PyG not installed")
def test_cgcnn_forward_pass() -> None:
    from torch_geometric.data import Data, Batch

    from tier1.ch10.cgcnn import CGCNN, GaussianBasis

    basis = GaussianBasis(0.0, 8.0, 16)
    n_atoms, n_edges = 4, 6
    Z = torch.tensor([6, 6, 8, 8], dtype=torch.long)
    edge_index = torch.tensor(
        [[0, 1, 1, 2, 2, 3], [1, 0, 2, 1, 3, 2]], dtype=torch.long
    )
    r = torch.full((n_edges,), 1.5)
    e = basis(r)
    data = Data(Z=Z, edge_index=edge_index, edge_attr=e)
    batch = Batch.from_data_list([data, data])
    model = CGCNN(atom_dim=16, edge_dim=16, n_conv=2, hidden_dim=32)
    out = model(batch)
    assert out.shape == (2,)


@pytest.mark.skipif(not torch_available(), reason="PyTorch / PyG not installed")
def test_cgcnn_training_reduces_loss() -> None:
    from torch_geometric.data import Data, Batch

    from tier1.ch10.cgcnn import CGCNN, GaussianBasis

    torch.manual_seed(0)
    basis = GaussianBasis(0.0, 8.0, 16)
    samples = []
    targets = []
    for i in range(8):
        Z = torch.tensor([6, 6, 8, 8], dtype=torch.long)
        edge_index = torch.tensor(
            [[0, 1, 1, 2, 2, 3], [1, 0, 2, 1, 3, 2]], dtype=torch.long
        )
        r = torch.full((6,), 1.0 + 0.1 * i)
        e = basis(r)
        y = torch.tensor([float(i)])
        samples.append(Data(Z=Z, edge_index=edge_index, edge_attr=e, y=y))
        targets.append(float(i))
    batch = Batch.from_data_list(samples)

    model = CGCNN(atom_dim=16, edge_dim=16, n_conv=2, hidden_dim=32)
    opt = torch.optim.Adam(model.parameters(), lr=5e-3)
    loss_fn = torch.nn.L1Loss()

    model.train()
    initial_loss = None
    final_loss = None
    for step in range(50):
        opt.zero_grad()
        pred = model(batch)
        loss = loss_fn(pred, batch.y.view(-1))
        if step == 0:
            initial_loss = float(loss)
        loss.backward()
        opt.step()
        final_loss = float(loss)
    assert final_loss is not None and initial_loss is not None
    assert final_loss < initial_loss, (
        f"Loss did not decrease: {initial_loss:.3f} -> {final_loss:.3f}"
    )
