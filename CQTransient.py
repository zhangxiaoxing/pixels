import numpy as np
import h5py

fstr = np.load("ctd.npz", allow_pickle=True)
features_per_su = fstr["features_per_su"].tolist()
reg_list = fstr["reg_list"].tolist()

# with h5py.File("su_trials_fr_6.hdf5", "w") as fw:
#     fw.create_dataset('count',data=len(features_per_su))
#     for suid in np.arange(0, len(features_per_su)):
#         S1_raw = features_per_su[suid]['S1_6']
#         S2_raw = features_per_su[suid]['S2_6']
#
#         S1 = np.array([np.mean(S1_raw[bin_idx:bin_idx + 4, :], axis=0) for bin_idx in np.arange(16, 40, 4)])
#         S2 = np.array([np.mean(S2_raw[bin_idx:bin_idx + 4, :], axis=0) for bin_idx in np.arange(16, 40, 4)])
#
#         grp = fw.create_group(suid.__str__())
#         grp.create_dataset('S1', data=S1)
#         grp.create_dataset('S2', data=S2)
#

def create_3s_dataset():
    with h5py.File("su_trials_fr_3.hdf5", "w") as fw:
        fw.create_dataset('count',data=len(features_per_su))
        for suid in np.arange(0, len(features_per_su)):
            S1_raw = features_per_su[suid]['S1_3']
            S2_raw = features_per_su[suid]['S2_3']

            S1 = np.array([np.mean(S1_raw[bin_idx:bin_idx + 4, :], axis=0) for bin_idx in np.arange(16, 28, 4)])
            S2 = np.array([np.mean(S2_raw[bin_idx:bin_idx + 4, :], axis=0) for bin_idx in np.arange(16, 28, 4)])

            grp = fw.create_group(suid.__str__())
            grp.create_dataset('S1', data=S1)
            grp.create_dataset('S2', data=S2)
