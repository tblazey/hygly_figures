import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Colors for each gene
colors = ['#006bb6', '#b6006b', '#00b64b', '#b64b00']

# Load in expression data
hk_non = pd.read_csv('Nonneurons_hk.csv')
hk_neu = pd.read_csv('Neurons_hk.csv')
hk_non_sc = pd.read_csv('Nonneurons_hk_superclass.csv', usecols=[1])

# Binarize gene expression matrices
hk_neu = hk_neu.astype({'HK1':bool, 'HK2':bool, 'HK3':bool, 'GCK':bool})
hk_non = hk_non.astype({'HK1':bool, 'HK2':bool, 'HK3':bool, 'GCK':bool})

# Get fraction of cells with some expression
hk_non_sums = hk_non.sum() / hk_non.shape[0]
hk_neu_sums = hk_neu.sum() / hk_neu.shape[0]

# Make neuron/non neuron plots
fig, ax = plt.subplots(1, 2, figsize=(14, 4), width_ratios=[0.75, 1])
bottom = np.array([0.0, 0.0])
ax[0].grid(linewidth=0.25, c='gray', which='major')
for idx, gene in enumerate(hk_non_sums.index):
    exp_frac = np.array([hk_neu_sums[gene], hk_non_sums[gene]])
    ax[0].bar(
        ['Neuronal', 'Non-Neuronal'],
        exp_frac,
        label=gene,
        bottom=bottom,
        color=colors[idx],
        edgecolor='black',
        linewidth=0.75,
        zorder=2
    )
    bottom += exp_frac
ax[0].set_xlabel('Cell Class', fontweight='bold', fontsize=12)
ax[0].set_ylabel('Fraction of Cells', fontweight='bold', fontsize=12)
ax[0].tick_params(axis='both', labelsize=10, width=0.5, colors='gray')
for axis in ['top','bottom','left','right']:
    ax[0].spines[axis].set_linewidth(0.5)
    ax[0].spines[axis].set_color('gray')
ax[0].text(-0.2, 1.05, 'C)', weight='bold', size=20, transform=ax[0].transAxes)

# Get fractions per cell type
hk_non['Class'] = hk_non_sc
cell_fracs = hk_non.groupby('Class').sum().div(hk_non['Class'].value_counts(), axis=0)
bottom = np.zeros(10)
for idx, gene in enumerate(hk_non_sums.index):
    ax[1].bar(
        cell_fracs.index,
        cell_fracs[gene],
        label=gene,
        bottom=bottom,
        color=colors[idx],
        edgecolor='black',
        linewidth=0.75,
        zorder=2
    )
    bottom += cell_fracs[gene]
ax[1].grid(linewidth=0.25, c='gray', which='major')
ax[1].set_xlabel('Cell Type', fontweight='bold', fontsize=12)
ax[1].set_ylabel('Fraction of Cells', fontweight='bold', fontsize=12)
ax[1].tick_params(axis='both', labelsize=10, width=0.5, colors='gray')
plt.setp(ax[1].xaxis.get_majorticklabels(), rotation=45, ha="right" )
for axis in ['top','bottom','left','right']:
    ax[1].spines[axis].set_linewidth(0.5)
    ax[1].spines[axis].set_color('gray')
ax[1].text(-0.2, 1.05, 'D)', weight='bold', size=20, transform=ax[1].transAxes)
ax[1].legend(
    bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0, fontsize=10
)

# Save combined figure
plt.subplots_adjust(wspace=0.4)
plt.savefig('figure_6cd.png', dpi=300, bbox_inches='tight')
plt.close('all')




