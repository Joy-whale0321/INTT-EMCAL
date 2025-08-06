import torch
from torch.utils.data import Dataset
import numpy as np
import joblib

from data_dphi import TrackCaloDataset as dphiDataset
from data_energy import TrackCaloDataset as energyDataset

from model_dphi import TrackCaloRegressor as dphiRegressor
from model_energy import TrackCaloRegressor as energyRegressor


class FusionDataset(Dataset):
    def __init__(self, list_file, pt_min=0.0, pt_max=10.0, device="cuda" if torch.cuda.is_available() else "cpu"):
        self.device = device

        # load dphi model
        model_dphi = dphiRegressor()
        model_dphi.load_state_dict(torch.load("/mnt/e/sphenix/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/version4/model_weight/best_model_pt_0.0_10.0_INTT_CaloIwoE.pt"))
        model_dphi.to(device)
        model_dphi.eval()

        # load energy model
        model_energy = energyRegressor()
        model_energy.load_state_dict(torch.load("/mnt/e/sphenix/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/version2/model_weight/best_model_pt_0.0_10.0_INTT_CaloIwoE.pt"))
        model_energy.to(device)
        model_energy.eval()

        scaler_energy = joblib.load("/mnt/e/sphenix/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/version2/model_weight/scaler_pt_0.0_10.0_INTT_CaloIwoE.pkl")

        # load datasets
        dataset_proxy = dphiDataset(list_file, pt_min=pt_min, pt_max=pt_max)

        dataset_trkcalo = energyDataset(list_file, pt_min=pt_min, pt_max=pt_max)
        dataset_trkcalo.X = torch.tensor(scaler_energy.transform(dataset_trkcalo.X.numpy()), dtype=torch.float32)

        assert len(dataset_proxy) == len(dataset_trkcalo)
        N = len(dataset_proxy)

        X_fusion = []
        Y_fusion = []

        with torch.no_grad():
            for i in range(N):
                x1 = dataset_proxy.X[i].unsqueeze(0).to(device)
                x2 = dataset_trkcalo.X[i].unsqueeze(0).to(device)
                
                dphi_i = dataset_proxy.X[i][-1]  # 最后一维
                calo_edep_i = dataset_trkcalo.X[i][-1]  # 最后一维
                pt_pred1_i = model_dphi(x1).cpu().numpy().item()
                pt_pred2_i = model_energy(x2).cpu().numpy().item()

                # how to conside pt bin (est? embed?)
                pt_est = 0.5 * (pt_pred1_i + pt_pred2_i)
                # embed
                pt_bin_edges = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]  # 11 个 edges → 10 个 bin
                pt_bin_onehot = [0] * 10  # 10 维 one-hot
                
                for j in range(len(pt_bin_edges)-1):
                    if pt_est >= pt_bin_edges[j] and pt_est < pt_bin_edges[j+1]:
                        pt_bin_onehot[j] = 1
                        break
                
                # 拼成 X_combined
                X_combined_i = [dphi_i, pt_pred1_i, calo_edep_i, pt_pred2_i] + pt_bin_onehot
                X_fusion.append(X_combined_i)

                # X_fusion.append([dphi_i , pt_pred1_i , calo_edep_i , pt_pred2_i, pt_est])

                y = dataset_proxy.Y[i]
                Y_fusion.append(y.item())

        self.X = torch.tensor(np.array(X_fusion), dtype=torch.float32)
        self.Y = torch.tensor(np.array(Y_fusion), dtype=torch.float32).view(-1, 1)

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return self.X[idx], self.Y[idx]
