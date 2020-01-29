import tensortools as tt
import numpy as np
import matplotlib.pyplot as plt

# Make synthetic dataset.


def nonneg_tca(X, R, prefix=''):

    # Fit CP tensor decomposition (two times).
    U = tt.ncp_bcd(X, rank=R, verbose=True)
    V = tt.ncp_bcd(X, rank=R, verbose=True)

    # Compare the low-dimensional factors from the two fits.
    # fig, ax, po = tt.plot_factors(U.factors)
    # tt.plot_factors(V.factors, fig=fig)
    # fig.suptitle("raw models")
    # fig.tight_layout()

    # Align the two fits and print a similarity score.
    sim = tt.kruskal_align(U.factors, V.factors, permute_U=True, permute_V=True)
    print(sim)

    # Plot the results again to see alignment.
    fig, ax, po = tt.plot_factors(U.factors, plots=["scatter", "scatter", "line"])
    tt.plot_factors(V.factors, plots=["scatter", "scatter", "line"], fig=fig)
    [x.set_xticks([11.5, 15.5, 27.5, 31.5]) for x in ax[:, 2]]

    ax[-1, 0].set_xlabel("SU #")
    ax[-1, 1].set_xlabel("Trial #")
    ax[-1, 2].set_xlabel("Time (s)")
    ax[-1, 2].set_xticklabels(["S", "+1", "T", "+1"])

    fig.suptitle("aligned models")
    fig.tight_layout()

    # Show plots.
    plt.show()
    fig.set_size_inches(40, 40)
    fig.set_dpi(300)
    fig.savefig(
        prefix+"nonneg_TCA_trial_" + str(X.shape[1]) + "_R" + str(R) + ".png",
        dpi=300,
        bbox_inches="tight",
    )
    return (U.obj, V.obj, sim)
