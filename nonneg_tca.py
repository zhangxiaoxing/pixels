import tensortools as tt
import numpy as np
import matplotlib.pyplot as plt

# Make synthetic dataset.

def nonneg_tca(X, R):
    
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
    fig, ax, po = tt.plot_factors(U.factors)
    tt.plot_factors(V.factors, fig=fig)
    fig.suptitle("aligned models")
    fig.tight_layout()
    
    # Show plots.
    plt.show()
    fig.set_size_inches(40,40)
    fig.set_dpi(300)
    fig.savefig('TCR_R'+str(R)+'.png',dpi=300,bbox_inches='tight')
    