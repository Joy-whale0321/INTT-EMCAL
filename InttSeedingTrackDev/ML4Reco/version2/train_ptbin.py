# === train.py ===
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split
from data import TrackCaloDataset
from model import TrackCaloRegressor
import os
import joblib

def train(list_file, pt_min=0.0, pt_max=2.0, batch_size=1024, epochs=300, lr=5e-5, val_ratio=0.3, device="cuda" if torch.cuda.is_available() else "cpu"):
    print(f"Training pt range: [{pt_min}, {pt_max}) GeV")
    print(f"Using device: {device}")

    dataset = TrackCaloDataset(list_file, pt_min=pt_min, pt_max=pt_max)
    train_size = int((1 - val_ratio) * len(dataset))
    val_size = len(dataset) - train_size
    train_set, val_set = random_split(dataset, [train_size, val_size])

    train_loader = DataLoader(train_set, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_set, batch_size=batch_size)

    model = TrackCaloRegressor().to(device)
    # optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-5)

    best_val_loss = float("inf")

    for epoch in range(epochs):
        model.train()
        train_loss = 0
        for xb, yb in train_loader:
            xb, yb = xb.to(device), yb.to(device)
            pred = model(xb)

            pt_reso = (yb - pred) / (yb)
            weights = (pt_reso.abs() < 0.2).float() * 2.0 + 1.0
            loss = ((pt_reso) ** 2 * weights).mean()
            # loss = (pt_reso ** 2).mean()

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * xb.size(0)
        train_loss /= len(train_loader.dataset)

        model.eval()
        val_loss = 0
        with torch.no_grad():
            for xb, yb in val_loader:
                xb, yb = xb.to(device), yb.to(device)
                pred = model(xb)
                pt_reso = (yb - pred) / (yb)
                weights = (pt_reso.abs() < 0.2).float() * 2.0 + 1.0
                loss = ((pt_reso) ** 2 * weights).mean()
                # loss = (pt_reso ** 2).mean()
                
                val_loss += loss.item() * xb.size(0)
        val_loss /= len(val_loader.dataset)

        print(f"Epoch {epoch+1:03d} | Train Loss: {train_loss:.4f} | Val Loss: {val_loss:.4f}")

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            torch.save(model.state_dict(), f"model_weight/best_model_pt_{pt_min:.1f}_{pt_max:.1f}_INTT_CaloIwoE.pt")
            print(f"✓ Saved best model (val loss = {val_loss:.4f})")

    print("✅ 训练完成。最优模型保存在 best_model_*.pt")

    joblib.dump(train_set.dataset.scaler, f"model_weight/scaler_pt_{pt_min:.1f}_{pt_max:.1f}_INTT_CaloIwoE.pkl")
    print("✅ 标准化器已保存")

if __name__ == "__main__":
    list_file = "../train500k.list"
    # pt_bins = [(0, 2), (2, 4), (4, 6), (6, 8), (8, 10)]
    pt_bins = [(0, 10)]
    for pt_min, pt_max in pt_bins:
        train(list_file, pt_min=pt_min, pt_max=pt_max)