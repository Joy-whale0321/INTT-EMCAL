import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split
from data import TrackCaloDataset
from model import TrackCaloRegressor
import os

def train(list_file, batch_size=1024, epochs=200, lr=5e-5, val_ratio=0.2, device="cuda" if torch.cuda.is_available() else "cpu"):
    print(f"Using device: {device}")

    dataset = TrackCaloDataset(list_file)
    train_size = int((1 - val_ratio) * len(dataset))
    val_size = len(dataset) - train_size
    train_set, val_set = random_split(dataset, [train_size, val_size])

    train_loader = DataLoader(train_set, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_set, batch_size=batch_size)

    model = TrackCaloRegressor().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    criterion = nn.MSELoss()

    best_val_loss = float("inf")

    for epoch in range(epochs):
        model.train()
        train_loss = 0
        for xb, yb in train_loader:
            xb, yb = xb.to(device), yb.to(device)
            pred = model(xb)

            pt_reso = (yb - pred) / yb
            loss = criterion(pt_reso, torch.zeros_like(pt_reso))
            # loss = criterion(pred, yb)

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

                pt_reso = (yb - pred) / yb
                loss = criterion(pt_reso, torch.zeros_like(pt_reso))
                # loss = criterion(pred, yb)

                val_loss += loss.item() * xb.size(0)
        val_loss /= len(val_loader.dataset)

        print(f"Epoch {epoch+1:03d} | Train Loss: {train_loss:.4f} | Val Loss: {val_loss:.4f}")

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            torch.save(model.state_dict(), "best_model.pt")
            print(f"✓ Saved best model (val loss = {val_loss:.4f})")

    torch.save(model.state_dict(), "final_model.pt")
    print("✅ 训练完成。最优模型保存在 best_model.pt, and final_model.pt")

if __name__ == "__main__":
    list_file = "../train50k.list"
    if not os.path.exists(list_file):
        print(f"❌ 找不到 {list_file}")
    else:
        train(list_file)
