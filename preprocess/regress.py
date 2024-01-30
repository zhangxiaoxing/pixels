import os, glob, sys
import pandas as pd
import numpy as np


def logistic_regress(path='.'):
    ks2dirs = glob.glob(r"imec?_ks2" + os.sep,root_dir=path)
    if len(ks2dirs) != 1:
        print("Unexpected imec_ks2 folder condition")
        sys.exit(101)
    metric_tbl = pd.read_csv(os.path.join(path,ks2dirs[0], "metrics.csv"))

    drift_metric = np.load(os.path.join(path,"drift_metrics.npy"))
    metric_tbl.insert(metric_tbl.shape[1], "drift_metric", drift_metric)

    curr_dir = os.path.dirname(__file__)
    qc_coeff = pd.read_csv(os.path.join(curr_dir, "qc_coeff.csv"))
    regress_val_arr = metric_tbl[
        [
            "firing_rate",
            "presence_ratio",
            "isi_viol",
            "amplitude_cutoff",
            "isolation_distance",
            "l_ratio",
            "d_prime",
            "nn_hit_rate",
            "nn_miss_rate",
            "silhouette_score",
            "max_drift",
            "cumulative_drift",
            "snr",
            "amplitude",
            "drift_metric",
        ]
    ].to_numpy()
    regress_coeff = (
        qc_coeff[["ALM", "Striatum", "Thalamus", "Midbrain", "Medulla"]]
        .iloc[0:15]
        .to_numpy()
    )
    qc_raw = (
        np.dot(regress_val_arr, regress_coeff)
        + qc_coeff[["ALM", "Striatum", "Thalamus", "Midbrain", "Medulla"]]
        .iloc[15]
        .to_numpy()
    )
    qc_logistic = 1 / (1 + np.exp(-1 * qc_raw))
    np.save(os.path.join(path,"logistic_regress.npy"), qc_logistic)
    return np.sum(qc_logistic[:, 0] > 0.5), qc_logistic.shape[0]


if __name__ == "__main__":
#     print(os.getcwd())
#     if os.path.isfile("logistic_regress.npy"):
#         print("data exist")
#         sys.exit(0)
# 
#     su_count, cluster_count = logistic_regress()
#     print(f"Number of SU={su_count} from {cluster_count} clusters")

    if "snakemake" in globals():
        outf=snakemake.output[0]
        print(f"Processing file: {outf}")
        su_count, cluster_count = logistic_regress(snakemake.wildcards[0])
        print(f"Number of SU={su_count} from {cluster_count} clusters")


