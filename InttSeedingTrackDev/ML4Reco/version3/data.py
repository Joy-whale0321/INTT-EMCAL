import numpy as np
import torch
from torch.utils.data import Dataset
import uproot
from collections import Counter
from numpy.linalg import norm
from sklearn.preprocessing import StandardScaler

def point_line_distance(point, line_point1, line_point2):
    point = np.array(point)
    line_point1 = np.array(line_point1)
    line_point2 = np.array(line_point2)
    line_vec = line_point2 - line_point1
    return norm(np.cross(point - line_point1, line_vec)) / norm(line_vec)

class TrackCaloDataset(Dataset):
    def __init__(self, list_file, tree_name="tree", scaler=None, pt_min=0.0, pt_max=10.0):
        data_X, data_Y = self.extract_features_from_rootlist(list_file, tree_name)

        data_X = np.array(data_X)
        data_Y = np.array(data_Y)
        mask = (data_Y >= pt_min) & (data_Y < pt_max)
        data_X = data_X[mask]
        data_Y = data_Y[mask]

        if scaler is None:
            scaler = StandardScaler()
            data_X = scaler.fit_transform(data_X)
        else:
            data_X = scaler.transform(data_X)

        self.X = torch.tensor(data_X, dtype=torch.float32)
        self.Y = torch.tensor(data_Y, dtype=torch.float32).view(-1, 1)
        self.scaler = scaler

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return self.X[idx], self.Y[idx]

    @staticmethod
    def extract_features_from_rootlist(list_file, tree_name="tree"):
        fail_012 = 0
        fail_34 = 0
        fail_56 = 0
        fail_calo = 0
        fail_truth = 0

        branches_to_load = [
            "trk_system", "trk_layer", "trk_X", "trk_Y", "trk_Z",
            "tower_system", "tower_X", "tower_Y", "tower_Z", "tower_edep",
            "PrimaryG4P_Pt"
        ]

        X_data = []
        Y_data = []

        with open(list_file, "r") as f:
            root_files = [line.strip() for line in f if line.strip()]

        for root_file in root_files:
            try:
                file = uproot.open(root_file)
                tree = file[tree_name]
                data = tree.arrays(branches_to_load, library="np")

                n_entries = len(data["trk_system"])

                for i in range(n_entries):
                    trk_layer = data["trk_layer"][i]
                    trk_x = data["trk_X"][i]
                    trk_y = data["trk_Y"][i]
                    trk_z = data["trk_Z"][i]
                    trk_hits = list(zip(trk_layer, trk_x, trk_y, trk_z))

                    clu_34 = [p for p in trk_hits if p[0] in (3, 4)]
                    clu_56 = [p for p in trk_hits if p[0] in (5, 6)]
                    if len(clu_34) != 1 or len(clu_56) != 1:
                        continue

                    p34 = np.array(clu_34[0][1:])
                    p56 = np.array(clu_56[0][1:])

                    track_point_layers = []
                    success = True
                    for layer_id in [0, 1, 2]:
                        layer_hits = [p for p in trk_hits if p[0] == layer_id]
                        if len(layer_hits) == 0:
                            success = False
                            break
                        dists = [point_line_distance(p[1:], p34, p56) for p in layer_hits]
                        min_idx = np.argmin(dists)
                        track_point_layers.append(layer_hits[min_idx][1:])

                    if not success:
                        continue

                    trk_feat = np.concatenate([
                        np.array(track_point_layers).flatten(),
                        p34,
                        p56
                    ])

                    tower_system = data["tower_system"][i]
                    tower_x = data["tower_X"][i]
                    tower_y = data["tower_Y"][i]
                    tower_z = data["tower_Z"][i]
                    tower_e = data["tower_edep"][i]

                    # 找到 system == 0 的数量作为 EMCal 数目
                    emcal_count = np.sum(tower_system == 0)

                    if emcal_count == 0:
                        fail_calo += 1
                        continue

                    tower_x = tower_x[:emcal_count]
                    tower_y = tower_y[:emcal_count]
                    tower_z = tower_z[:emcal_count]
                    tower_e = tower_e[:emcal_count]

                    tower_xyzE = np.stack([tower_x, tower_y, tower_z, tower_e], axis=1)

                    k = 9
                    valid_k = min(len(tower_xyzE), k)
                    if valid_k < k:
                        pad = np.zeros((k - valid_k, 4))
                        tower_xyzE = np.concatenate([tower_xyzE, pad], axis=0)
                    else:
                        idx = np.argsort(tower_xyzE[:, 3])[-k:]
                        tower_xyzE = tower_xyzE[idx]

                    tower_feat = tower_xyzE.flatten()
                    feat = np.concatenate([trk_feat, tower_feat, [valid_k]])
                    X_data.append(feat)

                    Truth_Pt = data["PrimaryG4P_Pt"][i]
                    if len(Truth_Pt) != 1:
                        continue
                    Y_data.append(Truth_Pt[0])

            except Exception as e:
                print(f"Error reading {root_file}: {e}")
                continue

        print(f"Total usable entries: {len(X_data)}")
        print(f"[Stats] Events failed cond_012: {fail_012}")
        print(f"[Stats] Events failed cond_3or4: {fail_34}")
        print(f"[Stats] Events failed cond_5or6: {fail_56}")

        return np.array(X_data), np.array(Y_data)
