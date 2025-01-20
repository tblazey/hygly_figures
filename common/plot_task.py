#!/usr/bin/python

#Load libraries
import os
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import nibabel as nib
import numpy as np
from scipy.ndimage import binary_dilation

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
base_min = 0
base_max = 100
diff_min = 0
diff_max = 75
gray_min = 3000
gray_max = 8000

#Get conditional means
eugly_mean = coef_data[:, :, :, 0]
hygly_mean = np.sum(coef_data[:, :, :, 0:2], axis=-1)
diff_mean = coef_data[..., 1]
eugly_mean[eugly_mean == 0] = np.nan
hygly_mean[hygly_mean == 0] = np.nan
diff_mean[diff_mean == 0] = np.nan
diff_thr = np.copy(diff_mean)
diff_thr[tfce_max < 0.95] = np.nan

# Dilate thresholded image so we can see boundary. Do it for each orthogonal slice
# so we aren't getting anything coming from another slice
diff_thr_mask = np.float32(~np.isnan(diff_thr))
diff_thr_dil_x = np.float32(binary_dilation(diff_thr_mask[44, :, :])) - diff_thr_mask[44, :, :]
diff_thr_dil_y = np.float32(binary_dilation(diff_thr_mask[:, 38, :])) - diff_thr_mask[:, 38, :]
diff_thr_dil_z = np.float32(binary_dilation(diff_thr_mask[:, :, 35])) - diff_thr_mask[:, :, 35]

# Mask out zeros
diff_thr_dil_x[diff_thr_dil_x == 0] = np.nan
diff_thr_dil_y[diff_thr_dil_y == 0] = np.nan
diff_thr_dil_z[diff_thr_dil_z == 0] = np.nan

#Sagittal hygly mean
ax[0, 0].matshow(
    np.rot90(brain_data[44, :, :]), cmap=gray_map, vmin=gray_min, vmax=gray_max
)
pos = ax[0, 0].matshow(
    np.rot90(hygly_mean[44, :, :]), cmap=autm_map, vmin=base_min, vmax=base_max
)
neg = ax[0, 0].matshow(
    np.rot90(hygly_mean[44, :, :]) * -1, cmap=cyan_map, vmin=base_min, vmax=base_max
)
ax[0, 0].axis('off')

#Sagittal eugly mean
ax[1, 0].matshow(
    np.rot90(brain_data[44, :, :]), cmap=gray_map, vmin=gray_min, vmax=gray_max
)
ax[1, 0].matshow(
    np.rot90(eugly_mean[44, :, :]), cmap=autm_map, vmin=base_min, vmax=base_max
)
ax[1, 0].matshow(
    np.rot90(eugly_mean[44, :, :]) * -1, cmap=cyan_map, vmin=base_min, vmax=base_max
)
ax[1, 0].axis('off')

#Sagittal diff mean
ax[2, 0].matshow(
    np.rot90(brain_data[44, :, :]), cmap=gray_map, vmin=gray_min, vmax=gray_max
)
ax[2, 0].matshow(
    np.rot90(diff_mean[44, ...]),
    cmap=autm_map,
    vmin=diff_min,
    vmax=diff_max,
    alpha=0.3
)
ax[2, 0].matshow(
    np.rot90(diff_mean[44, ...]) * -1,
    cmap=cyan_map,
    vmin=diff_min,
    vmax=diff_max,
    alpha=0.3
)
pos_diff = ax[2, 0].matshow(
    np.rot90(diff_thr[44, ...]), cmap=autm_map, vmin=diff_min, vmax=diff_max
)
neg_diff = ax[2, 0].matshow(
    np.rot90(diff_thr[44, ...]) * -1, cmap=cyan_map, vmin=diff_min, vmax=diff_max
)
ax[2, 0].imshow(np.rot90(diff_thr_dil_x), cmap='gray', vmax=1)
ax[2, 0].axis('off')

#Coronal eugly mean
ax[0, 1].matshow(
    np.rot90(brain_data[:, 38, :]), cmap=gray_map, vmin=gray_min, vmax=gray_max
)
ax[0, 1].matshow(
    np.rot90(hygly_mean[:, 38, :]), cmap=autm_map, vmin=base_min, vmax=base_max
)
ax[0, 1].matshow(
    np.rot90(hygly_mean[:, 38, :]) * -1, cmap=cyan_map, vmin=base_min, vmax=base_max
)
ax[0, 1].set_title('%s: Hygly.'%(task_name), fontweight='bold', fontsize=16, y=0.93)
ax[0, 1].axis('off')

#Coronal eugly mean
ax[1, 1].matshow(
    np.rot90(brain_data[:, 38, :]), cmap=gray_map, vmin=gray_min, vmax=gray_max
)
ax[1, 1].matshow(
    np.rot90(eugly_mean[:, 38, :]), cmap=autm_map, vmin=base_min, vmax=base_max
)
ax[1, 1].matshow(
    np.rot90(eugly_mean[:, 38, :]) * -1, cmap=cyan_map, vmin=base_min, vmax=base_max
)
ax[1, 1].set_title('Eugly.', fontweight='bold', fontsize=16, y=0.93)
ax[1, 1].axis('off')

#Coronal diff mean
ax[2, 1].matshow(
    np.rot90(brain_data[:, 38, :]), cmap=gray_map, vmin=gray_min, vmax=gray_max
)
ax[2, 1].matshow(
    np.rot90(diff_mean[:, 38, :]),
    cmap=autm_map,
    vmin=diff_min,
    vmax=diff_max,
    alpha=0.3
)
ax[2, 1].matshow(
    np.rot90(diff_mean[:, 38, :]) * -1,
    cmap=cyan_map,
    vmin=diff_min,
    vmax=diff_max,
    alpha=0.3
)
ax[2, 1].matshow(
    np.rot90(diff_thr[:, 38, :]), cmap=autm_map, vmin=diff_min, vmax=diff_max
)
ax[2, 1].matshow(
    np.rot90(diff_thr[:, 38, :]) * -1, cmap=cyan_map, vmin=diff_min, vmax=diff_max,
)
ax[2, 1].set_title('Hygly. - Eugly.', fontweight='bold', fontsize=16, y=0.93)
ax[2, 1].imshow(np.rot90(diff_thr_dil_y), cmap='gray', vmax=1)
ax[2, 1].axis('off')

#Axial eugly mean
ax[0, 2].matshow(
    np.rot90(brain_data[:, :, 35]), cmap=gray_map, vmin=gray_min, vmax=gray_max
)
ax[0, 2].matshow(
    np.rot90(hygly_mean[:, :, 35]), cmap=autm_map, vmin=base_min, vmax=base_max
)
ax[0, 2].matshow(
    np.rot90(hygly_mean[:, :, 35]) * -1, cmap=cyan_map, vmin=base_min, vmax=base_max
)
ax[0, 2].axis('off')

#Axial eugly mean
ax[1, 2].matshow(
    np.rot90(brain_data[:, :, 35]), cmap=gray_map, vmin=gray_min, vmax=gray_max
)
ax[1, 2].matshow(
    np.rot90(eugly_mean[:, :, 35]), cmap=autm_map, vmin=base_min, vmax=base_max
)
ax[1, 2].matshow(
    np.rot90(eugly_mean[:, :, 35]) * -1, cmap=cyan_map, vmin=base_min, vmax=base_max
)
ax[1, 2].axis('off')

#Axial diff mean
ax[2, 2].matshow(
    np.rot90(brain_data[:, :, 35]), cmap=gray_map, vmin=gray_min, vmax=gray_max
)
ax[2, 2].matshow(
    np.rot90(diff_mean[..., 35]),
    cmap=autm_map,
    vmin=diff_min,
    vmax=diff_max,
    alpha=0.3
)
ax[2, 2].matshow(
    np.rot90(diff_mean[..., 35]) * -1,
    cmap=cyan_map,
    vmin=diff_min,
    vmax=diff_max,
    alpha=0.3
)
ax[2, 2].matshow(
    np.rot90(diff_thr[..., 35]), cmap=autm_map, vmin=diff_min, vmax=diff_max,
)
ax[2, 2].matshow(
    np.rot90(diff_thr[..., 35]) * -1, cmap=cyan_map, vmin=diff_min, vmax=diff_max
)
ax[2, 2].imshow(np.rot90(diff_thr_dil_z), cmap='gray', vmax=1)
ax[2, 2].axis('off')

#Add colorbars
pos_bar_ax = fig.add_axes([0.90, 0.625, 0.03, 0.2])
neg_bar_ax = fig.add_axes([0.9, 0.4, 0.03, 0.2])
neg_bar_ax.invert_yaxis()
fig.colorbar(pos, cax=pos_bar_ax)
fig.colorbar(neg, cax=neg_bar_ax)

# Add colorbars for difference
pos_bar_diff_ax = fig.add_axes([0.90, 0.265, 0.03, 0.1])
neg_bar_diff_ax = fig.add_axes([0.9, 0.14, 0.03, 0.1])
neg_bar_diff_ax.invert_yaxis()
fig.colorbar(pos_diff, cax=pos_bar_diff_ax)
fig.colorbar(neg_diff, cax=neg_bar_diff_ax)

#Add colorbar labels
pos_bar_ax.annotate('Task > Rest', [2.5, 0.35], fontweight='bold',
                    xycoords='axes fraction', rotation=90)
neg_bar_ax.annotate('Task < Rest', [2.5, 0.35], fontweight='bold',
                    xycoords='axes fraction', rotation=90) 
neg_bar_diff_ax.annotate('COPE', [-0.2, -0.25], fontweight='bold', xycoords='axes fraction')         
                  

#Add subplot labels
ax[0, 0].annotate("A)", [0, 1], fontweight='bold', xycoords='axes fraction', fontsize=16)
ax[1, 0].annotate("B)", [0, 1], fontweight='bold', xycoords='axes fraction', fontsize=16)                           
ax[2, 0].annotate("C)", [0, 1], fontweight='bold', xycoords='axes fraction', fontsize=16)

plt.subplots_adjust(hspace=0.2, wspace=0.01)
plt.savefig(f'{figure}.tiff', dpi=300, bbox_inches='tight')
plt.close('all')
