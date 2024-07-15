#!/usr/bin/python

#Load libraries
import os
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import nibabel as nib
import numpy as np

# Setup
brain_path = os.path.join(
    os.environ['FSLDIR'], 'data/standard/MNI152_T1_2mm_brain.nii.gz'
)
if sys.argv[1] == 'enc':
    task = 'enc'
    task_name = 'Encoding'
    figure = 'figure_s7'
elif sys.argv[1] == 'ret':
    task = 'ret'
    task_name = 'Retrieval'
    figure = 'figure_s8'
else:
    sys.exit()

#Load in coefficients
coef_hdr = nib.load(f'all_feat1_{task}_glm_pe.nii.gz')
coef_data = coef_hdr.get_fdata()

#Load in all p-values
tfce_1 = nib.load(f'all_feat1_{task}_tfce_corrp_tstat1.nii.gz')
tfce_2 = nib.load(f'all_feat1_{task}_tfce_corrp_tstat2.nii.gz')
tfce_3 = nib.load(f'all_feat2_{task}_tfce_corrp_tstat1.nii.gz')
tfce_4 = nib.load(f'all_feat2_{task}_tfce_corrp_tstat2.nii.gz')
tfce_1_data = tfce_1.get_fdata()
tfce_2_data = tfce_2.get_fdata()
tfce_3_data = tfce_3.get_fdata()
tfce_4_data = tfce_4.get_fdata()

#Make max p-value image
tfce_max = np.max(np.stack((tfce_1_data, tfce_2_data, tfce_3_data, tfce_4_data)), axis=0)

#Load in fsl brain image
brain_hdr = nib.load(brain_path)
brain_data = brain_hdr.get_fdata()

#Define colormaps
autm_map = plt.get_cmap('autumn').copy()
autm_map.set_under(alpha=0)
cyan_data = autm_map(np.arange(255))[:, [2, 1, 0, 3]]
cyan_map = ListedColormap(cyan_data)
cyan_map.set_under(alpha=0)
gray_map = plt.get_cmap('gray').copy()
gray_map.set_under('white', alpha=0)

#Create figure (slices are: 44, 38, 35)
fig, ax = plt.subplots(3, 3, figsize=(10, 10))
base_min = 15
base_max = 100
diff_min = 10
diff_max = 75
gray_min = 3000
gray_max = 8000

#Get conditional means
eugly_mean = coef_data[:, :, :, 0]

#Sagittal eugly mean
ax[0, 0].matshow(
    np.rot90(brain_data[44, :, :]), cmap=gray_map, vmin=gray_min, vmax=gray_max
)
pos = ax[0, 0].matshow(
    np.rot90(eugly_mean[44, :, :]), cmap=autm_map, vmin=base_min, vmax=base_max
)
neg = ax[0, 0].matshow(
    np.rot90(eugly_mean[44, :, :]) * -1, cmap=cyan_map, vmin=base_min, vmax=base_max
)
ax[0, 0].axis('off')

#Sagittal diff mean
ax[1, 0].matshow(
    np.rot90(brain_data[44, :, :]), cmap=gray_map, vmin=gray_min, vmax=gray_max
)
ax[1, 0].matshow(
    np.rot90(coef_data[44, :, :, 1]), cmap=autm_map, vmin=diff_min, vmax=diff_max
)
ax[1, 0].matshow(
    np.rot90(coef_data[44, :, :, 1]) * -1, cmap=cyan_map, vmin=diff_min, vmax=diff_max
)
ax[1, 0].axis('off')

#Sagittal p-value
ax[2, 0].matshow(
    np.rot90(brain_data[44, :, :]), cmap=gray_map, vmin=gray_min, vmax=gray_max
)
p_val = ax[2, 0].matshow(np.rot90(tfce_max[44, :, :]), cmap=autm_map, vmin=0.95, vmax=1)
ax[2, 0].axis('off')

#Coronal eugly mean
ax[0, 1].matshow(
    np.rot90(brain_data[:, 38, :]), cmap=gray_map, vmin=gray_min, vmax=gray_max
)
ax[0, 1].matshow(
    np.rot90(eugly_mean[:, 38, :]), cmap=autm_map, vmin=base_min, vmax=base_max
)
ax[0, 1].matshow(
    np.rot90(eugly_mean[:, 38, :]) * -1, cmap=cyan_map, vmin=base_min, vmax=base_max
)
ax[0, 1].set_title('%s: Eugly.'%(task_name), fontweight='bold', fontsize=16, y=0.93)
ax[0, 1].axis('off')

#Coronal diff mean
ax[1, 1].matshow(
    np.rot90(brain_data[:, 38, :]), cmap=gray_map, vmin=gray_min, vmax=gray_max
)
ax[1, 1].matshow(
    np.rot90(coef_data[:, 38, :, 1]), cmap=autm_map, vmin=diff_min, vmax=diff_max
)
ax[1, 1].matshow(
    np.rot90(coef_data[:, 38, :, 1]) * -1, cmap=cyan_map, vmin=diff_min, vmax=diff_max
)
ax[1, 1].set_title('Hygly. - Eugly.', fontweight='bold', fontsize=16, y=0.93)
ax[1, 1].axis('off')

#Coonal p-value
ax[2, 1].matshow(
    np.rot90(brain_data[:, 38, :]), cmap=gray_map, vmin=gray_min, vmax=gray_max
)
ax[2, 1].matshow(np.rot90(tfce_max[:, 38, :]), cmap=autm_map, vmin=0.95, vmax=1)
ax[2, 1].set_title('TFCE', fontweight='bold', fontsize=16, y=0.93)
ax[2, 1].axis('off')

#Axial eugly mean
ax[0, 2].matshow(
    np.rot90(brain_data[:, :, 35]), cmap=gray_map, vmin=gray_min, vmax=gray_max
)
ax[0, 2].matshow(
    np.rot90(eugly_mean[:, :, 35]), cmap=autm_map, vmin=base_min, vmax=base_max
)
ax[0, 2].matshow(
    np.rot90(eugly_mean[:, :, 35]) * -1, cmap=cyan_map, vmin=base_min, vmax=base_max
)
ax[0, 2].axis('off')

#Axial diff mean
ax[1, 2].matshow(
    np.rot90(brain_data[:, :, 35]), cmap=gray_map, vmin=gray_min, vmax=gray_max
)
ax[1, 2].matshow(
    np.rot90(coef_data[:, :, 35, 1]), cmap=autm_map, vmin=diff_min, vmax=diff_max
)
ax[1, 2].matshow(
    np.rot90(coef_data[:, :, 35, 1]) * -1, cmap=cyan_map, vmin=diff_min, vmax=diff_max
)
ax[1, 2].axis('off')

#Axial p-value
ax[2, 2].matshow(
    np.rot90(brain_data[:, :, 35]), cmap=gray_map, vmin=gray_min, vmax=gray_max
)
ax[2, 2].matshow(np.rot90(tfce_max[:, :, 35]), cmap=autm_map, vmin=0.95, vmax=1)
ax[2, 2].axis('off')

#Add colorbars
pos_bar_ax = fig.add_axes([0.90, 0.625, 0.03, 0.2])
neg_bar_ax = fig.add_axes([0.9, 0.4, 0.03, 0.2])
p_bar_ax = fig.add_axes([0.9, 0.1, 0.03, 0.2])
neg_bar_ax.invert_yaxis()
fig.colorbar(pos, cax=pos_bar_ax)
fig.colorbar(neg, cax=neg_bar_ax)
fig.colorbar(p_val, cax=p_bar_ax)
p_bar_ax.axes.set_yticks([0.95, 1.0], labels=[0.05, 1E-5])

#Add colorbar labels
pos_bar_ax.annotate('Task > Rest', [2.5, 0.35], fontweight='bold',
                    xycoords='axes fraction', rotation=90)
neg_bar_ax.annotate('Task < Rest', [2.5, 0.35], fontweight='bold',
                    xycoords='axes fraction', rotation=90)     
p_bar_ax.annotate('p-value', [2.5, 0.35], fontweight='bold',
                  xycoords='axes fraction', rotation=90)   
                  

#Add subplot labels
ax[0, 0].annotate("A)", [0, 1], fontweight='bold', xycoords='axes fraction', fontsize=16)
ax[1, 0].annotate("B)", [0, 1], fontweight='bold', xycoords='axes fraction', fontsize=16)                           
ax[2, 0].annotate("C)", [0, 1], fontweight='bold', xycoords='axes fraction', fontsize=16)

plt.subplots_adjust(hspace=0.2, wspace=0.01)
plt.savefig(f'{figure}.tiff', dpi=300, bbox_inches='tight')
plt.close('all')
