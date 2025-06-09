import torch
import ROOT
import numpy as np
import joblib
from model import TrackCaloRegressor
from data import TrackCaloDataset

device = "cuda" if torch.cuda.is_available() else "cpu"
pt_bins = [(0, 2), (2, 4), (4, 6), (6, 8), (8, 10), (0, 10)]
# pt_bins = [(0, 2)]

out_file = ROOT.TFile("outputFile/pt_relative_error_Si_CaloTower.root", "RECREATE")

for pt_min, pt_max in pt_bins:
    print(f"\n=== Testing model for pt in [{pt_min}, {pt_max}) GeV ===")

    # 加载 scaler 和数据（注意要传 pt 范围）
    # scaler = joblib.load(f"model_weight/scaler_pt_{pt_min:.1f}_{pt_max:.1f}.pkl")
    scaler = joblib.load(f"model_weight/scaler_pt_0.0_10.0.pkl")
    dataset = TrackCaloDataset("../test500k.list", scaler=scaler, pt_min=pt_min, pt_max=pt_max)
    loader = torch.utils.data.DataLoader(dataset, batch_size=256)

    # 加载模型
    model = TrackCaloRegressor()
    # model.load_state_dict(torch.load(f"model_weight/best_model_pt_{pt_min:.1f}_{pt_max:.1f}.pt"))
    model.load_state_dict(torch.load(f"model_weight/best_model_pt_0.0_10.0.pt"))
    model.to(device)
    model.eval()

    # 直方图
    name = f"h_relerr_{pt_min}_{pt_max}"
    title = f"RelErr for {pt_min} < pt < {pt_max};(pred - truth)/truth;Counts"
    hist = ROOT.TH1D(name, title, 200, -1, 1)

    # 推理并填充
    with torch.no_grad():
        for x, y in loader:
            x = x.to(device)
            y = y.to(device)
            pred = model(x)

            y_np = y.cpu().numpy().flatten()
            pred_np = pred.cpu().numpy().flatten()
            rel_err = (pred_np - y_np) / (y_np + 1e-6)

            for rei in rel_err:
                hist.Fill(rei)

    hist.Fit("gaus", "", "", -0.1, 0.1)
    hist.Write()

out_file.Close()
print("✅ 每个 pt 区间的相对误差已保存为 pt_relative_error_bybin.root")
