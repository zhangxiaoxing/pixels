import sys, os, glob
import numpy as np
import pandas as pd
from scipy.stats import poisson


def drift_metric(SpkID, SpkTS, trials, suids):
    driftwin = np.stack((trials[:, 1] - 5 * 30000, trials[:, 1]), axis=1)
    drift_metric_list = []
    for oneid in suids:
        oneTS = SpkTS[SpkID == oneid]
        inwin = (oneTS[:, None] > driftwin[:, 0]) & (oneTS[:, None] < driftwin[:, 1])
        counts = np.sum(inwin, axis=0)
        mm = np.mean(counts)
        (lb, ub) = (poisson.ppf(0.025, mm), poisson.ppf(0.975, mm))
        outlier_count = np.sum(np.logical_or(counts < lb, counts > ub))
        drift_metric_list.append(outlier_count / counts.shape[0])
    drift_array = np.array(drift_metric_list, dtype=np.float32)
    return drift_array


def load_data():
    trials = np.load("sync_trials.npy")
    if trials.shape[0] < 80:
        print("Unexpected trials format")
        sys.exit(102)
    ks2dirs = glob.glob(r"imec?_ks2" + os.sep)
    if len(ks2dirs) != 1:
        print("Unexpected imecx_ks2 folder condition")
        sys.exit(101)
    SpkTS = np.load(os.path.join(ks2dirs[0], "spike_times.npy"))
    SpkID = np.load(os.path.join(ks2dirs[0], "spike_clusters.npy"))

    metric_tbl = pd.read_csv(os.path.join(ks2dirs[0], "metrics.csv"))
    suids = metric_tbl["cluster_id"].to_numpy()

    return SpkTS, SpkID, trials, suids


if __name__ == "__main__":  # Debuging entry. Import module when possible.
    print(os.getcwd())
    if os.path.exists("drift_metric.npy"):
        print("data exist")
        sys.exit(0)
    if not os.path.exists("sync_trials.npy"):
        print("Missing data")
        sys.exit(103)
    (SpkTS, SpkID, trials, suids) = load_data()
    drift_metric_arr = drift_metric(SpkID, SpkTS, trials, suids)
    np.save("drift_metric.npy", drift_metric_arr)
