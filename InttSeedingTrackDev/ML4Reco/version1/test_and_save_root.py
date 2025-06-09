import torch
import ROOT
from model import TrackCaloRegressor
from data import TrackCaloDataset
import numpy as np

# ==== 加载数据集 ====
dataset = TrackCaloDataset("../test50k.list")
loader = torch.utils.data.DataLoader(dataset, batch_size=128)
device = "cuda" if torch.cuda.is_available() else "cpu"

# ==== 加载模型 ====
model = TrackCaloRegressor()
model.load_state_dict(torch.load("best_model.pt"))
# model.load_state_dict(torch.load("final_model.pt"))
model.to(device)
model.eval()

# ==== 创建 ROOT 直方图 ====
h1 = ROOT.TH1D("h_relative_error", "Relative Error: (pred - truth)/truth;Rel Error;Counts", 400, -2, 2)

# ==== 推理并填充 ====
with torch.no_grad():
    for x, y in loader:
        x = x.to(device)
        y = y.to(device)
        pred = model(x)
        # 计算 relative error
        rel_err = ((pred - y) / y).cpu().numpy().flatten()
        for val in rel_err:
            h1.Fill(val)

# ==== 高斯拟合并保存到 ROOT 文件 ====
fit_result = h1.Fit("gaus", "", "", -0.08, 0.15)  # "S" 返回 TF1 对象
fit_func = h1.GetFunction("gaus")

out_file = ROOT.TFile("pt_relative_error.root", "RECREATE")
h1.Write()
if fit_func:
    fit_func.Write()  # 保存拟合曲线
out_file.Close()

print("✅ 相对误差直方图及高斯拟合已写入 pt_relative_error.root")
