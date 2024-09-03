import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd

#Change default math font
mpl.rcParams['mathtext.default'] = 'regular'

# Load in atlas image
atlas_path = './wmparc_raw_comb.nii.gz'
atlas_hdr = nib.load(atlas_path)
atlas_img = atlas_hdr.get_fdata().flatten()

# Find unique regions
rois = np.unique(atlas_img)[1::]

# Load in regional data 
reg_df =  pd.read_excel('../common/data_values_file.xlsx', 'Allen.Gene')

# Make empty image for showing regional expression
hk1_img = np.zeros_like(atlas_img)
hk2_img = np.zeros_like(atlas_img)
ratio_img = np.zeros_like(atlas_img)

# Make expression images
for idx, roi in enumerate(reg_df.id):
   if np.isnan(reg_df['HK1'][idx]) == False:
      hk1_img[atlas_img == roi] = reg_df['HK1'][idx]
      hk2_img[atlas_img == roi] =  reg_df['HK2'][idx]
      ratio_img[atlas_img == roi] =  reg_df['Ratio'][idx]

# Reshape images and set zeros to nan
hk1_img = hk1_img.reshape(atlas_hdr.shape)
hk1_img[hk1_img == 0] = np.nan
hk2_img = hk2_img.reshape(atlas_hdr.shape)
hk2_img[hk2_img == 0] = np.nan
ratio_img = ratio_img.reshape(atlas_hdr.shape)
ratio_img[ratio_img == 0] = np.nan

# Make expression image figure
fig, ax = plt.subplots(1, 3, figsize=(16, 8), gridspec_kw={'width_ratios':[1, 1, 1]})

# Add HK1 image
hk = ax[0].matshow(np.rot90(hk1_img[25:170, 25:205, 96]), cmap='plasma')
ax[0].set_title('HK1', fontweight='bold', fontsize=28)
ax[0].text(-0.2, 1.2, 'A)', weight='bold', size=36, transform=ax[0].transAxes)
ax[0].axis('off')

# Add HK2 image
ax[1].matshow(np.rot90(hk2_img[25:170, 25:205, 96]), cmap='plasma')
ax[1].axis('off')
ax[1].set_title('HK2', fontweight='bold', fontsize=28)

# Colorbar for HK images
cb_ax_1 = fig.add_axes([0.26, 0.125, 0.225, 0.04])
c_bar_1 = fig.colorbar(hk, cax=cb_ax_1, orientation='horizontal')
c_bar_1.set_label('Norm. Exp.', fontsize=18, labelpad=15, fontweight='bold')

# Add ratio image
ratio = ax[2].matshow(np.rot90(ratio_img[25:170, 25:205, 96]), cmap='plasma')
ax[2].axis('off')
ax[2].set_title('HK1 / HK2', fontweight='bold', fontsize=28)

# Colorbar for ratio
cb_ax = fig.add_axes([0.66, 0.125, 0.225, 0.04])
c_bar = fig.colorbar(ratio, cax=cb_ax, orientation='horizontal')
c_bar.set_label('HK 1 / HK 2', fontsize=18, labelpad=15, fontweight='bold')
plt.subplots_adjust(wspace=0.15)

plt.savefig('figure_6a.png', dpi=300, bbox_inches='tight')
plt.close('all')

# Make scatter plot figure
fig2, ax2 = plt.subplots(1, 1, figsize=(4, 4))
ax2.grid(linewidth=0.25, c='gray')
for name, group in reg_df.groupby('Group'):
    ax2.plot(group.Ratio, group.Delta, marker='o', ms=4, ls='', label=name)

# Format plot
ax2.legend(
    loc="upper center",
    bbox_to_anchor=(0.5, 1.1),
    borderaxespad=0,
    fontsize=8,
    ncol=np.unique(reg_df.Group).shape[0],
    columnspacing=0.7,
    handlelength=0.2,
    frameon=False
)
ax2.set_xlabel('HK1 / HK2 Ratio', fontweight='bold', fontsize=10)
ax2.set_ylabel(
    r'$\Delta$ [$^{18}$F]FDG (SUVR)', fontweight='bold', fontsize=10, labelpad=0
)
ax2.tick_params(axis='both', labelsize=8, width=0.5, colors='gray')
for axis in ['top','bottom','left','right']:
    ax2.spines[axis].set_linewidth(0.5)
    ax2.spines[axis].set_color('gray')
ax2.text(-0.3, 1, 'B)', weight='bold', size=18, transform=ax2.transAxes)
cor = np.corrcoef(reg_df['Ratio'], reg_df['Delta'])[0, 1]
print(np.corrcoef(reg_df['Ratio'], reg_df['Basal'])[0, 1])
ax2.text(0.7, 0.9, f'r = {cor:.2f}', weight='bold', size=8.5, transform=ax2.transAxes)

plt.savefig('figure_6b.png', dpi=300, bbox_inches='tight')
plt.close('all')




