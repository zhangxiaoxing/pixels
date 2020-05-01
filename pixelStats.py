import numpy as np

def baselineVector(one_su_trial_FR):
    base = one_su_trial_FR[2:10, :].flatten()
    if np.std(base):
        return (np.mean(base), np.std(base))
    else:
        base = one_su_trial_FR[-10:-1, :].flatten()
        if np.std(base):
            return (np.mean(base), np.std(base))
        else:
            print('flat base')
            base = one_su_trial_FR.flatten()
            return (0, np.std(base))