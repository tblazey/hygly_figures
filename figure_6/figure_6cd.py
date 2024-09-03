import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Colors for each gene
colors = ['#006bb6', '#b6006b', '#00b64b', '#b64b00']

# Load in Siletti expression data
gm_hk_non = pd.read_csv('Nonneurons_hk.csv')
gm_hk_neu = pd.read_csv('Neurons_hk.csv')
gm_hk_non_sc = pd.read_csv('Nonneurons_hk_superclass.csv', usecols=[1])

# Binarize gene expression matrices
gm_hk_neu = gm_hk_neu.astype({'HK1':bool, 'HK2':bool, 'HK3':bool, 'GCK':bool})
gm_hk_non = gm_hk_non.astype({'HK1':bool, 'HK2':bool, 'HK3':bool, 'GCK':bool})

# Load in Seeker expression data
wm_hk = pd.read_csv('./white_matter_hk_data.csv')

# Get fraction of cells/class that express each gene in each cell type
wm_cell_fracs = wm_hk.drop(columns=['tissue', 'class']).groupby('cell_type').sum().div(wm_hk['cell_type'].value_counts(), axis=0)
wm_class_fracs = wm_hk.drop(columns=['tissue', 'cell_type']).groupby('class').sum().div(wm_hk['class'].value_counts(), axis=0)
wm_class_fracs = wm_class_fracs.drop(columns='HK3')

# Load in cell type labels
wm_cell_lbls = pd.read_csv('white_matter_hk_types.csv')
wm_cell_fracs = pd.merge(wm_cell_fracs, wm_cell_lbls, on='cell_type', how='left')
wm_cell_fracs = wm_cell_fracs.drop(columns='HK3')
wm_hk2_frac = np.sum(wm_hk.HK2) / wm_hk.shape[0]
wm_hk2_class_frac = np.sum(wm_hk['class'][wm_hk.HK2] == 'Non-Neuronal') / wm_hk['class'][wm_hk.HK2].shape[0]
print(wm_hk2_frac * 100)
print(wm_hk2_class_frac * 100)
print(wm_class_fracs * 100)

# Remove neuronal cells
wm_non_neuron_mask = np.logical_and.reduce((
    wm_cell_fracs.cell_type != 'glutamatergic neuron',
    wm_cell_fracs.cell_type != 'GABAergic neuron',
    wm_cell_fracs.cell_type != 'cerebellar granule cell',
    wm_cell_fracs.cell_type != 'neuron'
))
wm_cell_fracs = wm_cell_fracs.loc[wm_non_neuron_mask]

# Get fraction of cells with some expression
gm_hk_non_sums = gm_hk_non.sum() / gm_hk_non.shape[0]
gm_hk_neu_sums = gm_hk_neu.sum() / gm_hk_neu.shape[0]
gm_class_fracs = pd.DataFrame({'Neuronal':gm_hk_neu_sums, 'Non-Neuronal':gm_hk_non_sums}).T
gm_hk2_frac = (np.sum(gm_hk_neu.HK2) + np.sum(gm_hk_non.HK2)) / (gm_hk_neu.shape[0] + gm_hk_non.shape[0])
gm_hk2_class_frac = np.sum(gm_hk_non.HK2) / (np.sum(gm_hk_neu.HK2) + np.sum(gm_hk_non.HK2))
print(gm_hk2_frac * 100)
print(gm_hk2_class_frac * 100)
print(gm_class_fracs * 100)

# Get fractions per cell type
gm_hk_non['Class'] = gm_hk_non_sc
gm_cell_fracs = gm_hk_non.groupby('Class').sum().div(gm_hk_non['Class'].value_counts(), axis=0)

# Make neuron/non neuron plots
fig, ax = plt.subplots(1, 2, figsize=(14, 4))

# Make Siletti plot
bottom = np.zeros(gm_cell_fracs.shape[0])
for idx, gene in enumerate(gm_hk_non_sums.index):
    ax[0].bar(
        gm_cell_fracs.index,
        gm_cell_fracs[gene],
        label=gene,
        bottom=bottom,
        color=colors[idx],
        edgecolor='black',
        linewidth=0.75,
        zorder=2
    )
    bottom += gm_cell_fracs[gene]
ax[0].grid(linewidth=0.25, c='gray', which='major')
ax[0].set_xlabel('Cell Type', fontweight='bold', fontsize=12, labelpad=-10)
ax[0].set_ylabel('Fraction of Cells', fontweight='bold', fontsize=12)
ax[0].set_ylim([0, 1])
ax[0].tick_params(axis='both', labelsize=10, width=0.5, colors='gray')
plt.setp(ax[0].xaxis.get_majorticklabels(), rotation=45, ha="right" )
for axis in ['top','bottom','left','right']:
    ax[0].spines[axis].set_linewidth(0.5)
    ax[0].spines[axis].set_color('gray')
ax[0].text(-0.2, 1.05, 'C)', weight='bold', size=18.6, transform=ax[0].transAxes)
ax[0].legend(
    bbox_to_anchor=(0.8, 0.8), loc="center left", borderaxespad=0, fontsize=10
)
ax[0].set_title('Siletti et al. 2023', fontweight='bold', fontsize=14, y=1.05)

# Seeker plot
colors.pop(2)
bottom = np.zeros(wm_cell_fracs.shape[0])
for idx, gene in enumerate(wm_class_fracs.columns):
    ax[1].bar(
        wm_cell_fracs.cell_type_labels,
        wm_cell_fracs[gene],
        label=gene,
        bottom=bottom,
        color=colors[idx],
        edgecolor='black',
        linewidth=0.75,
        zorder=2
    )
    bottom += wm_cell_fracs[gene]
ax[1].grid(linewidth=0.25, c='gray', which='major')
ax[1].set_xlabel('Cell Type', fontweight='bold', fontsize=12, labelpad=-30)
ax[1].set_ylabel('Fraction of Cells', fontweight='bold', fontsize=12)
ax[1].set_ylim([0, 1])
ax[1].tick_params(axis='both', labelsize=10, width=0.5, colors='gray')
plt.setp(ax[1].xaxis.get_majorticklabels(), rotation=45, ha="right" )
for axis in ['top','bottom','left','right']:
    ax[1].spines[axis].set_linewidth(0.5)
    ax[1].spines[axis].set_color('gray')
ax[1].text(-0.2, 1.05, 'D)', weight='bold', size=18.6, transform=ax[1].transAxes)
ax[1].legend(
    bbox_to_anchor=(0.8, 0.8), loc="center left", borderaxespad=0, fontsize=10
)
ax[1].set_title('Seeker et al. 2023', fontweight='bold', fontsize=14, y=1.05)

# Save combined figure
plt.subplots_adjust(wspace=0.3)
plt.savefig('figure_6cd.png', dpi=300, bbox_inches='tight')
plt.close('all')




