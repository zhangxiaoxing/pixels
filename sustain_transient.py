import numpy as np

# %% load file

fstr = np.load('GLM_stats.npz')
reg_arr = fstr['reg_arr']
all_sess_arr = fstr['all_sess_arr']

sel_wedge = [np.count_nonzero(all_sess_arr[x, :]) for x in np.arange(0, 6)]
sum_sel = np.sum(sel_wedge)

sum_nonmod = np.sum(np.count_nonzero(all_sess_arr[12, :]))

any_sel = np.any(all_sess_arr[:6, :], axis=0)

np.count_nonzero(any_sel & all_sess_arr[12, :].astype(np.bool))
