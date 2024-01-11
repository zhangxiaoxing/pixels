import glob, re, os, warnings, json, sys
import numpy as np
from scipy.stats import linregress

sps = 30000
# os.chdir(os.path.expanduser("~/npdata_out"))
# all_prb_trls = glob.glob("**/sync_trials.npy", recursive=True)
#
# sess = set()
# for onetrl in all_prb_trls:
#     match = re.search(r"catgt_.{5,40}?(?=/.*)", onetrl)
#     if match:
#         sess.add(match.group())
#     else:
#         warnings.warn("Error processing data folder " + onetrl)

# for one_sess in sess:
#     prb_sel = list(one_sess in trials for trials in all_prb_trls)
#     probes = np.sort(np.compress(prb_sel, all_prb_trls))

# sess_path = re.findall(r"^.*?_g\d" + os.sep, ref_trls_file)[0]
sess_path = snakemake.wildcards[0]
print("Processing " + sess_path)

# existing file preferably handled by snakemake
# if os.path.isfile(os.path.join(sess_path, "sess_spk_t_id.npy")):
#     print("Session spk time and id exist")
#     sys.exit(0)

probes = sorted(
    glob.glob(os.path.join(sess_path, "**/sync_trials.npy"), recursive=True)
)
ref_trls_file = probes[0]

ref_trls = np.load(ref_trls_file)
ref_idx = next(re.finditer(r"(?<=imec)\d", ref_trls_file))[0]
ref_spk_file = os.path.join(
    ref_trls_file.replace("sync_trials.npy", "imec" + ref_idx + "_ks2"),
    "spike_times.npy",
)
if not os.path.isfile(ref_spk_file):
    print("Missing file " + ref_spk_file)
    sys.exit(1)
ref_spk_t = np.load(ref_spk_file) / sps
ref_spk_id = np.load(ref_spk_file.replace("spike_times.npy", "spike_templates.npy"))

aligned_spk_t = np.column_stack((np.zeros(ref_spk_t.shape), ref_spk_t, ref_spk_id))
raw_spk_t = np.column_stack((np.ones(ref_spk_t.shape), ref_spk_t, ref_spk_id))

trls_match = np.ones((ref_trls.shape[0], 1))
trls_idx = np.full((ref_trls.shape[0], 1), ref_idx, dtype=float)
sess_trials = np.hstack((trls_idx, trls_match, ref_trls))

if len(probes) > 1:
    regress_dict = {}
    goodmatch = True
    for probe in probes[1:]:
        curr_idx = next(re.finditer(r"(?<=imec)\d", probe))[0]
        curr_trls_file = probe
        curr_trls = np.load(curr_trls_file)
        curr_spk_file = os.path.join(
            curr_trls_file.replace("sync_trials.npy", "imec" + curr_idx + "_ks2"),
            "spike_times.npy",
        )
        curr_spk_ts = np.load(curr_spk_file)
        curr_spk_id = np.load(
            curr_spk_file.replace("spike_times.npy", "spike_templates.npy")
        )

        probe_spk_t = curr_spk_ts / sps
        raw_spk_t = np.concatenate(
            (
                raw_spk_t,
                np.column_stack(
                    (
                        np.full(probe_spk_t.shape, curr_idx, dtype=int),
                        probe_spk_t,
                        curr_spk_id + int(curr_idx) * 10000,
                    )
                ),
            ),
            axis=0,
        )

        trls_match = np.zeros((curr_trls.shape[0], 1))
        trls_idx = np.full((curr_trls.shape[0], 1), curr_idx, dtype=float)
        sess_trials = np.vstack(
            (sess_trials, np.hstack((trls_idx, trls_match, curr_trls)))
        )

        if curr_trls.shape == ref_trls.shape:
            onset_delta = curr_trls[:, 0] - ref_trls[:, 0]
            qtrs = np.percentile(onset_delta, [25, 50, 75])
            iqr = qtrs[2] - qtrs[0]
            outlier = np.abs(onset_delta - qtrs[1]) > (iqr * 1.5)
            if np.any(outlier):
                print(sess_path)
                print(f"outliers in trials num={np.sum(outlier)}")
            if np.sum(outlier) > 5:
                goodmatch = False
                print(f"Bad match")
                break
            lr = linregress(
                curr_trls[np.logical_not(outlier), 0],
                ref_trls[np.logical_not(outlier), 0],
            )
            regress_dict[ref_idx + "_" + curr_idx] = lr
            if np.square(lr.rvalue) < 0.9999999:
                print(f"Low quality fit r={lr.rvalue}")
                goodmatch = False
                break
            aligned_t = (curr_spk_ts * lr.slope + lr.intercept) / sps
            aligned_spk_t = np.concatenate(
                (
                    aligned_spk_t,
                    np.column_stack(
                        (
                            np.zeros(aligned_t.shape),
                            aligned_t,
                            curr_spk_id + int(curr_idx) * 10000,
                        )
                    ),
                ),
                axis=0,
            )
            prb_trl_match = np.array(np.logical_not(outlier), dtype=float)
            sess_trials[sess_trials[:, 0] == float(curr_idx), 1] = prb_trl_match

    if goodmatch:
        np.save(os.path.join(sess_path, "sess_spk_t_id.npy"), aligned_spk_t)
    with open(os.path.join(sess_path, f"sess_time_regress.json"), "w") as f:
        json.dump(regress_dict, f)

if len(probes) == 1 or not goodmatch:
    np.save(os.path.join(sess_path, "sess_spk_t_id.npy"), raw_spk_t)

np.save(os.path.join(sess_path, "sess_trls.npy"), sess_trials)

sys.exit(0)
