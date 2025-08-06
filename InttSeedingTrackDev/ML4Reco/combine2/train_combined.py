import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split

from model_combined import FusionRegressor
from data_combined import FusionDataset

def train_fusion(list_file, pt_min=0.0, pt_max=10.0, batch_size=1024, epochs=300, device="cuda" if torch.cuda.is_available() else "cpu"):
    print(f"=== Training FusionRegressor: pt range [{pt_min}, {pt_max}) GeV ===")
    print(f"Using device: {device}")

    dataset = FusionDataset(list_file, pt_min=pt_min, pt_max=pt_max)

    # split train / val
    train_size = int(0.7 * len(dataset))
    val_size = len(dataset) - train_size
    train_set, val_set = random_split(dataset, [train_size, val_size])

    train_loader = DataLoader(train_set, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_set, batch_size=batch_size)

    model = FusionRegressor().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=5e-5, weight_decay=1e-5)

    best_val_loss = float("inf")

    for epoch in range(epochs):
        # === train ===
        model.train()
        train_loss = 0
        for xb, yb in train_loader:
            xb, yb = xb.to(device), yb.to(device)
            pred = model(xb)
            pt_reso = (yb - pred) / (yb)
            weights = (pt_reso.abs() < 0.2).float() * 2.0 + 1.0
            loss = ((pt_reso) ** 2 * weights).mean()

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * xb.size(0)
        train_loss /= len(train_loader.dataset)

        # === val ===
        model.eval()
        val_loss = 0
        with torch.no_grad():
            for xb, yb in val_loader:
                xb, yb = xb.to(device), yb.to(device)
                pred = model(xb)
                pt_reso = (yb - pred) / (yb)
                weights = (pt_reso.abs() < 0.2).float() * 2.0 + 1.0
                loss = ((pt_reso) ** 2 * weights).mean()

                val_loss += loss.item() * xb.size(0)
        val_loss /= len(val_loader.dataset)

        print(f"Epoch {epoch+1:03d} | Train Loss: {train_loss:.4f} | Val Loss: {val_loss:.4f}")

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            torch.save(model.state_dict(), "model_weight/best_model_combined.pt")
            print(f"✓ Saved best combined model (val loss = {val_loss:.4f})")

    print("✅ Training complete. Best model saved as best_model_combined.pt")


if __name__ == "__main__":
    list_file = "../train50k.list"
    pt_bins = [(0, 10.0)]
    for pt_min, pt_max in pt_bins:
        train_fusion(list_file, pt_min=pt_min, pt_max=pt_max)
