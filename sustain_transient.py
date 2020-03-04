import numpy as np
import matplotlib.pyplot as plt

# %% load file

fstr = np.load('GLM_stats.npz')
reg_arr = fstr['reg_arr']
all_sess_arr = fstr['all_sess_arr']

sel_wedge = [np.count_nonzero(all_sess_arr[x, :]) for x in np.arange(0, 6)]
sum_sel = np.sum(sel_wedge)

sum_nonmod = np.sum(np.count_nonzero(all_sess_arr[12, :]))

any_sel = np.any(all_sess_arr[:6, :], axis=0)

# so some of the su are considered non-modulated since the activity are comparable to BS, however S1 S2 trial are
# separated in statistics
_awkward = np.count_nonzero(any_sel & all_sess_arr[12, :].astype(np.bool))

frac = [sum_sel, all_sess_arr.shape[1] - sum_sel - sum_nonmod, sum_nonmod]
explode = (0.2, 0, 0)
labels = ('Sample selective', 'Non-selective modulation', 'Unmodulated')

(fh, axes) = plt.subplots(1, 2, figsize=(12, 6), dpi=50)
axes[0].pie(frac, explode=explode, labels=labels, autopct='%1.1f%%', shadow=True)
axes[0].axis('equal')

labels = ('Sample selective', 'Non-selective modulation', 'Unmodulated')
explode = (0, 0, 0, 0, 0, 0.2)
labels = ('only sample', 'only early delay', 'only late delay', 'SP & ED', 'ED & LD', 'SP & ED & LD')

axes[1].pie(sel_wedge, explode=explode, labels=labels, autopct=lambda p: '{:.1f}%'.format(p * 0.212), radius=0.5,
            shadow=True)
axes[1].axis('equal')
axes[0].set_xlim((-1.25, 1.25))
axes[1].set_xlim((-1.25, 1.25))
plt.show()
