import torch
import ROOT
import numpy as np

from model_combined import FusionRegressor
from data_combined import FusionDataset

device = "cuda" if torch.cuda.is_available() else "cpu"

pt_bins = [(0.0, 10.0)]

out_file = ROOT.TFile("outputFile/pt_relative_error_combined.root", "RECREATE")

for pt_min, pt_max in pt_bins:
    print(f"\n=== Testing CombinedRegressor for pt in [{pt_min}, {pt_max}) GeV ===")

    dataset = FusionDataset("../test500k.list", pt_min=pt_min, pt_max=pt_max)
    loader = torch.utils.data.DataLoader(dataset, batch_size=512)

    model = FusionRegressor()
    model.load_state_dict(torch.load("model_weight/best_model_combined.pt"))
    model.to(device)
    model.eval()

    # 2D hist: pt vs relative error
    nbins_2dy = (pt_max - pt_min) * 10
    hist2d = ROOT.TH2D("h2_pt_vs_relerr_combined", "pt vs RelErr Combined;truth pt;(pred - truth)/truth", int(nbins_2dy), float(pt_min), float(pt_max), 250, -2., 2.)

    # 1D hist: relative error
    hist1d = ROOT.TH1D("h_relerr_combined", "RelErr Combined;(pred - truth)/truth;Counts", 200, -1, 1)

    with torch.no_grad():
        for xb, yb in loader:
            xb, yb = xb.to(device), yb.to(device)
            pred = model(xb)

            y_np = yb.cpu().numpy().flatten()
            pred_np = pred.cpu().numpy().flatten()
            # Apply post-correction
            for i in range(len(pred_np)):
                reco_pt = pred_np[i]

                # Apply correction for reco_pt >= 8.8
                if reco_pt >= 8.8:
                    correction_factor = 0.02 + 0.08 * (reco_pt - 8.8)
                    pred_np[i] = reco_pt * (1.0 + correction_factor)

            rel_err = (pred_np - y_np) / (y_np)

            for yi, rei in zip(y_np, rel_err):
                hist1d.Fill(rei)
                hist2d.Fill(yi, rei)

    # # 高斯拟合
    # fit_result = hist1d.Fit("gaus", "S", -0.1, 0.1)
    # gaus_func = hist1d.GetFunction("gaus")
    # mean = gaus_func.GetParameter(1)
    # sigma = gaus_func.GetParameter(2)
    # print(f"Gaussian fit result: mean = {mean:.4f}, sigma = {sigma:.4f}")

    # c1 = ROOT.TCanvas("c1", "RelErr Combined", 800, 600)
    # hist1d.Draw()
    # gaus_func.SetLineColor(ROOT.kRed)
    # gaus_func.SetLineWidth(2)
    # gaus_func.Draw("same")

    # c1.SaveAs("outputFile/h_relerr_combined_fit.png")

    hist1d.Write()
    hist2d.Write()

out_file.Close()
print("✅ CombinedRegressor test result saved to pt_relative_error_combined.root")
