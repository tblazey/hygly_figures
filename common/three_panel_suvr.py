import sys
import fastkde
import matplotlib as mpl
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
from PIL import Image
from scipy.interpolate import interpn

#Change default math font
mpl.rcParams['mathtext.default'] = 'regular'

# Constants
    
if sys.argv[1] == 'rogi':
    figure = 'figure_4'
    mod_label = 'rOGI'
    img_min = 0.8
    img_max = 1.4
    diff_max = 0.5
elif sys.argv[1] == 'roef':
    figure = 'figure_s3'
    mod_label = 'rOEF'
    img_min = 0.8
    img_max = 1.2
    diff_max = 0.25
elif sys.argv[1] == 'fdg':
    figure = 'figure_s4'
    mod_label = r'[$^{18}F$]FDG'
    img_min = 0.5
    img_max = 1.5
    diff_max = 0.25
elif sys.argv[1] == 'co':
    figure = 'figure_s5'
    mod_label = r'[$^{15}O$]CO'
    img_min = 0.5
    img_max = 3.0
    diff_max = 0.25
else:
    sys.exit()
mod = f'{sys.argv[1]}_suvr'
cb_label = 'SUVR'

if sys.argv[1] == 'rogi' or sys.argv[1] == 'roef':
    plt_title = fr'{mod_label}: Hyper'
else:
    plt_title = fr'{mod_label} Uptake: Hyper'

img_mid = (img_max - img_min) / 2.0 + img_min
x = 48
y = 54
z = 38
dpi = 300.0
diff_min = 0
diff_mid = (diff_max - diff_min) / 2.0 + diff_min

#Load in image headers
eugly = nib.load(f'{mod}_basal.nii.gz')
hygly = nib.load(f'{mod}_hypergly.nii.gz')
diff = nib.load(f'{mod}_hypergly_coef.nii.gz')
fdr = nib.load(f'{mod}_hypergly_logp_fdr_05.nii.gz')
mpr = nib.load(f'../common/template_2mm.nii.gz')
mask = nib.load(f'../common/template_2mm_mask.nii.gz')
roi = nib.load(f'../common/MNI152_wmparc_comb_2mm_on_icbm_dwm.nii.gz')

#Load in image data
eugly_data = np.squeeze(eugly.get_fdata())
hygly_data = np.squeeze(hygly.get_fdata())
fdr_data = np.squeeze(fdr.get_fdata())
diff_data = np.squeeze(diff.get_fdata())
mpr_data = np.float32(np.squeeze(mpr.get_fdata()))
mask_data = np.logical_or(np.squeeze(mask.get_fdata()) == 0.0, eugly_data == 0.0)
roi_data = np.float32(np.squeeze(roi.get_fdata()))

# Get correlation values before masking
mask_invert = np.logical_not(mask_data)
base_corr = np.corrcoef(
    eugly_data[mask_invert].flatten(), diff_data[mask_invert].flatten())
print(base_corr[0, 1])

#Mask outside of fdr threshold
diff_data[np.abs(fdr_data) < 1.3] = np.nan

#Apply brain mask
eugly_data[mask_data] = np.nan
hygly_data[mask_data] = np.nan
diff_data[mask_data] = np.nan
mpr_data[mask_data] = np.nan
roi_data[mask_data] = np.nan
roi_data[roi_data == 0] = np.nan

#Make colormap for average images
img_map = plt.get_cmap("plasma")
img_map.set_bad("white", alpha=0)
img_map.set_under(img_map(0))
img_map.set_over(img_map(255))

#Make positive colormap
pos_map = plt.get_cmap('autumn')
pos_map.set_bad("white", alpha=0)
pos_map.set_under("white", alpha=0)
pos_map.set_over(pos_map(255), alpha=1)

#Make negative colormap
neg_lut = np.flipud(pos_map(np.arange(256))[:, 0:3]) * -1 + 1
neg_map = mpl.colors.ListedColormap(neg_lut)
neg_map.set_bad("white", alpha=0)
neg_map.set_under("white", alpha=0)
neg_map.set_over(neg_map(255), alpha=1)

#Make colorbar for negative colormap
disp_lut = pos_map(np.arange(256))[:, 0:3] * -1 + 1
disp_map = mpl.colors.ListedColormap(disp_lut)

#Make structural colormap
mpr_map = plt.get_cmap("gray")
mpr_map.set_bad("white", alpha=0)

#Make colormap for roi
roi_map = mpl.colors.ListedColormap([(0.75, 0, 0.43)])

#Make figure object
fig1 = plt.figure(figsize=[1500.0 / dpi, 600.0 / dpi])

#Make first grid
sp1 = gridspec.GridSpec(ncols=5, nrows=1, figure=fig1,
                        width_ratios=[218.0/182.0, 1, 1, 0.05, 0.05], wspace=0.05)

#Add subplots to grid
ax1 = fig1.add_subplot(sp1[0, 0])
ax2 = fig1.add_subplot(sp1[0, 1])
ax3 = fig1.add_subplot(sp1[0, 2])
ax4 = fig1.add_subplot(sp1[0, 3])

#Make first row
ax1.imshow(np.rot90(mpr_data[x, :, :]), cmap=mpr_map)
ax1.imshow(np.rot90(hygly_data[x, :, :]), cmap=img_map, vmin=img_min, vmax=img_max)
ax2.imshow(np.rot90(mpr_data[:, :, z]), cmap=mpr_map)
ax2.imshow(np.rot90(hygly_data[:, :, z]), cmap=img_map, vmin=img_min, vmax=img_max)
ax3.imshow(np.rot90(mpr_data[:, y, :]), cmap=mpr_map)
ax3.imshow(np.rot90(hygly_data[:, y, :]), cmap=img_map, vmin=img_min, vmax=img_max)

#Make first row
pax1 = ax1.get_position()
pax2 = ax2.get_position()
pax3 = ax3.get_position()
ax1.set_position([0.1, 0.125, pax1.width, pax1.height])
ax2.set_position([0.38, 0.08, pax2.width, pax2.height])
ax3.set_position([0.615, 0.12, pax3.width, pax3.height])

#Add title
ax1.text(0.1, 0.79, 'A)', weight='bold', size=8, transform=fig1.transFigure)
ax2.set_title(plt_title, weight='bold', size=8)

#Add colorbar
cb = mpl.colorbar.ColorbarBase(ax4, cmap=img_map, ticks=[0, 0.5, 1])
cb.ax.set_position([0.875, 0.33, 0.0225, 0.375])
cb.set_ticklabels(['%.1f'%(i) for i in [img_min, img_mid, img_max]])
for tick in cb.ax.yaxis.get_major_ticks():
    tick.label2.set_weight('bold')
    tick.label2.set_size(4.5)
cb.ax.set_xlabel(cb_label, size=5.5, weight='bold', labelpad=5)

#Turn off the axis labels
ax1.set_axis_off()
ax2.set_axis_off()
ax3.set_axis_off()

#Save figure
plt.savefig(f'{figure}a.png', dpi=dpi, bbox_inches='tight')
plt.close()

#Make second figure
fig2 = plt.figure(figsize=[1500.0 / dpi, 600.0 / dpi])

#Make second grid
sp2 = gridspec.GridSpec(ncols=5, nrows=1, figure=fig2,
                        width_ratios=[218.0/182.0, 1, 1, 0.05, 0.05], wspace=0.05)

#Add subplots to grid
ax6 = fig2.add_subplot(sp2[0, 0])
ax7 = fig2.add_subplot(sp2[0, 1])
ax8 = fig2.add_subplot(sp2[0, 2])
ax9 = fig2.add_subplot(sp2[0, 3])

#Make second row
ax6.imshow(np.rot90(mpr_data[x, :, :]), cmap=mpr_map)
ax6.imshow(np.rot90(eugly_data[x, :, :]), cmap=img_map, vmin=img_min, vmax=img_max)
ax7.imshow(np.rot90(mpr_data[:, :, z]), cmap=mpr_map)
ax7.imshow(np.rot90(eugly_data[:, :, z]), cmap=img_map, vmin=img_min, vmax=img_max)
ax8.imshow(np.rot90(mpr_data[:, y, :]), cmap=mpr_map)
ax8.imshow(np.rot90(eugly_data[:, y, :]), cmap=img_map, vmin=img_min, vmax=img_max)

#Make second row
pax6 = ax6.get_position()
pax7 = ax7.get_position()
pax8 = ax8.get_position()
ax6.set_position([0.1, 0.125, pax6.width, pax6.height])
ax7.set_position([0.38, 0.08, pax7.width, pax7.height])
ax8.set_position([0.615, 0.12, pax8.width, pax8.height])

#Add title
ax6.text(0.1, 0.79, 'B)', weight='bold', size=8, transform=fig2.transFigure)
ax7.set_title('Eugly.', weight='bold', size=8)

#Add colorbar
cb2 = mpl.colorbar.ColorbarBase(ax9, cmap=img_map, ticks=[0, 0.5, 1])
cb2.ax.set_position([0.875, 0.33, 0.0225, 0.375])
cb2.set_ticklabels(['%.1f'%(i) for i in [img_min, img_mid, img_max]])
for tick in cb2.ax.yaxis.get_major_ticks():
    tick.label2.set_weight('bold')
    tick.label2.set_size(4.5)
cb2.ax.set_xlabel(cb_label, size=5.5, weight='bold', labelpad=5)

#Turn off the axis labels
ax6.set_axis_off()
ax7.set_axis_off()
ax8.set_axis_off()

#Save figure
plt.savefig(f'{figure}b.png', dpi=dpi, bbox_inches='tight')
plt.close()

#Make last grid
fig3 = plt.figure(figsize=[1500.0 / dpi, 600.0 / dpi])
sp3 = gridspec.GridSpec(ncols=5, nrows=1, figure=fig3,
                        width_ratios=[218.0/182.0, 1, 1, 0.05, 0.05], wspace=0.05)

#Add subplots to grid
ax9 = fig3.add_subplot(sp3[0, 0])
ax10 = fig3.add_subplot(sp3[0, 1])
ax11= fig3.add_subplot(sp3[0, 2])
ax12 = fig3.add_subplot(sp3[0, 3])
ax13 = fig3.add_subplot(sp3[0, 4])

#Make first row
ax9.imshow(np.rot90(mpr_data[x, :, :]), cmap=mpr_map)
ax9.imshow(np.rot90(diff_data[x, :, :]), cmap=pos_map, vmin=diff_min, vmax=diff_max)
ax9.imshow(np.rot90(-diff_data[x, :, :]), cmap=neg_map, vmin=diff_min, vmax=diff_max)
ax10.imshow(np.rot90(mpr_data[:, :, z]), cmap=mpr_map)
ax10.imshow(np.rot90(diff_data[:, :, z]), cmap=pos_map, vmin=diff_min, vmax=diff_max)
ax10.imshow(np.rot90(-diff_data[:, :, z]), cmap=neg_map, vmin=diff_min, vmax=diff_max)
ax11.imshow(np.rot90(mpr_data[:, y, :]), cmap=mpr_map)
ax11.imshow(np.rot90(diff_data[:, y, :]), cmap=pos_map, vmin=diff_min, vmax=diff_max)
ax11.imshow(np.rot90(diff_data[:, y, :]), cmap=pos_map, vmin=diff_min, vmax=diff_max)
ax11.imshow(np.rot90(-diff_data[:, y, :]), cmap=neg_map, vmin=diff_min, vmax=diff_max)

#Move first row
pax1 = ax1.get_position()
pax2 = ax2.get_position()
pax3 = ax3.get_position()
ax9.set_position([0.1, 0.125, pax1.width, pax1.height])
ax10.set_position([0.38, 0.08, pax2.width, pax2.height])
ax11.set_position([0.615, 0.12, pax3.width, pax3.height])

#Add title
ax9.text(0.1, 0.79, 'C)', weight='bold', size=8, transform=fig3.transFigure)
ax10.set_title('Hyper. - Eugly.', weight='bold', size=8)

#Add positive colorbar
cb2 = mpl.colorbar.ColorbarBase(ax12, cmap=pos_map, ticks=[0, 0.5, 1])
cb2.ax.set_position([0.875, 0.57, 0.0225, 0.175])
cb2.set_ticklabels([f'{i:.2f}' for i in [diff_min, diff_mid, diff_max]])
for tick in cb2.ax.yaxis.get_major_ticks():
    tick.label2.set_weight('bold')
    tick.label2.set_size(4.5)

#Add negative colorbar
cb3 = mpl.colorbar.ColorbarBase(ax13, cmap=disp_map, ticks=[0, 0.5, 1])
cb3.ax.set_position([0.875, 0.33, 0.0225, 0.175])
cb3.set_ticklabels([f'{i:.2f}' for i in [-diff_max, -diff_mid, -diff_min]])
for tick in cb3.ax.yaxis.get_major_ticks():
    tick.label2.set_weight('bold')
    tick.label2.set_size(4.5)
cb3.ax.set_xlabel(cb_label, size=5.5, weight='bold', labelpad=3)

#Turn off the axis labels
ax9.set_axis_off()
ax10.set_axis_off()
ax11.set_axis_off()

#Save figure
plt.savefig(f'{figure}c.png', dpi=dpi, bbox_inches='tight')
plt.close()

#Read in each image
hygly_img = Image.open(f'{figure}a.png')
eug_img = Image.open(f'{figure}b.png')
diff_img = Image.open(f'{figure}c.png')

#Crop voxelwise images
coords = (0, 0, diff_img.size[0], diff_img.size[1] * 0.9)
hygly_crop = hygly_img.crop(coords)
eugly_crop = eug_img.crop(coords)
diff_crop = diff_img.crop(coords)

#Make new image
output = Image.new("RGBA",
                   (diff_crop.size[0], int(diff_crop.size[1] * 3)),
                   color=(255, 255, 255))

#Paste images to output
output.paste(hygly_crop, (0, 0))
output.paste(eugly_crop, (0, eugly_crop.size[1]))
output.paste(diff_crop, (0, eugly_crop.size[1] * 2))

#Save image
output.save(f'{figure}.tiff')