
# -*- coding: utf-8 -*-



def bars():
    # load and align 6s 3s delay selectivity data

    # plot transient bars
    # 3s transient from
    norm = Normalize(vmin=0, vmax=0.6)
    counts_sum = []
    (fig, ax) = plt.subplots(1, 1, figsize=(7.5, 7), dpi=300)
    for row in range(7):
        (counts, sums) = subgroup_equiv(3, row)
        counts_sum.extend(counts)
        ax.scatter(range(counts.shape[0]), np.ones(counts.shape[0]) * row, s=counts / 3, c=counts / sums, cmap='jet',
                   norm=norm)
    ax.set_xlim((-1, 7))
    ax.set_ylim((-1, 7))

    ax.set_xticks(np.arange(0, 7))
    ax.set_yticks(np.arange(0, 7))

    ax.set_yticklabels(('sustained', 'transient', 'switched', 'unclassified', 'only during sample',
                        'nonselective modulation', 'non-modulated'), rotation=45, va='top', ha='right')
    ax.set_xticklabels(('sustained', 'transient', 'switched', 'unclassified', 'only during sample',
                        'nonselective modulation', 'non-modulated'), rotation=45, va='top', ha='right')
    ax.set_xlabel('coding in 6s delay trials')
    ax.set_ylabel('coding in 3s delay trials')
    sm = plt.cm.ScalarMappable(cmap='jet', norm=norm)
    sm._A = []
    plt.colorbar(sm, ticks=[0, 0.5], format="%.1f")
    fig.savefig('sus_trans_3vs6.png', bbox_inches='tight', pad_inches=2)
    plt.show()
    sumcounts = np.sum(counts_sum)
    print(f"sum counts {sumcounts}")

    counts_sum = []
    (fig, ax) = plt.subplots(1, 1, figsize=(7.5, 7), dpi=300)
    for row in range(7):
        (counts, sums) = subgroup_equiv('early3in6', row)
        counts_sum.extend(counts)
        ax.scatter(range(counts.shape[0]), np.ones(counts.shape[0]) * row, s=counts / 3, c=counts / sums, cmap='jet',
                   norm=norm)
    ax.set_xlim((-1, 7))
    ax.set_ylim((-1, 7))

    ax.set_xticks(np.arange(0, 7))
    ax.set_yticks(np.arange(0, 7))

    ax.set_yticklabels(('sustained', 'transient', 'switched', 'unclassified', 'only during sample',
                        'nonselective modulation', 'non-modulated'), rotation=45, va='top', ha='right')
    ax.set_xticklabels(('sustained', 'transient', 'switched', 'unclassified', 'only during sample',
                        'nonselective modulation', 'non-modulated'), rotation=45, va='top', ha='right')
    ax.set_xlabel('coding in 6s delay early half')
    ax.set_ylabel('coding in 3s delay trials')
    sm = plt.cm.ScalarMappable(cmap='jet', norm=norm)
    sm._A = []
    plt.colorbar(sm, ticks=[0, 0.5], format="%.1f")
    fig.savefig('sus_trans_3vse3.png', bbox_inches='tight', pad_inches=2)
    plt.show()
    sumcounts = np.sum(counts_sum)
    print(f"sum counts {sumcounts}")

    counts_sum = []
    (fig, ax) = plt.subplots(1, 1, figsize=(7.5, 7), dpi=300)
    for row in range(7):
        (counts, sums) = subgroup_equiv('early_late', row)
        counts_sum.extend(counts)
        ax.scatter(range(counts.shape[0]), np.ones(counts.shape[0]) * row, s=counts / 3, c=counts / sums, cmap='jet',
                   norm=norm)
    ax.set_xlim((-1, 7))
    ax.set_ylim((-1, 7))

    ax.set_xticks(np.arange(0, 7))
    ax.set_yticks(np.arange(0, 7))

    ax.set_yticklabels(('sustained', 'transient', 'switched', 'unclassified', 'only during sample',
                        'nonselective modulation', 'non-modulated'), rotation=45, va='top', ha='right')
    ax.set_xticklabels(('sustained', 'transient', 'switched', 'unclassified', 'only during sample',
                        'nonselective modulation', 'non-modulated'), rotation=45, va='top', ha='right')
    ax.set_xlabel('coding in 6s delay early half')
    ax.set_ylabel('coding in 6s delay late half')

    sm = plt.cm.ScalarMappable(cmap='jet', norm=norm)
    sm._A = []
    plt.colorbar(sm, ticks=[0, 0.5], format="%.1f")
    fig.savefig('sus_trans_early_vs_late.png', bbox_inches='tight', pad_inches=2)
    plt.show()
    sumcounts = np.sum(counts_sum)
    print(f"sum counts {sumcounts}")