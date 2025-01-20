import os
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
from scipy.ndimage import binary_dilation

#Change default math font
mpl.rcParams['mathtext.default'] = 'regular'

# Constants
suff = ''
if sys.argv[1] == 'cbf':
    cb_label = '\\frac{mL}{hg \cdot min}'
    mod = 'cbf'
    mod_label = 'CBF'
    figure = 'figure_s6'
    img_min = 10
    img_max = 100 
    diff_min = 0
    diff_max = 10
    sx =  [0, 110]
    sy = [-25, 25]
    sy_ticks = [-20, -10, 0, 10, 20]
    wy = [5, 50]
    plt_title = rf'CBF: Eugly.'
    reg_mod = 'cbf'
elif sys.argv[1][0:6] == 'cmrglc':
    if sys.argv[1] == 'cmrglc':
        figure = 'figure_2'
        suff = ''
    else:
        figure = 'figure_s1'
        suff = '_lc'
    mod = 'cmrglc'
    cb_label =  '\\frac{\mu Mol}{hg \cdot min}'
    mod_label = 'CMRglc'
    img_min = 10    
    img_max = 50
    sx =  [10, 50]
    sy = [-10, 20]
    diff_min = 0
    diff_max = 15
    sy_ticks = [-10, -5, 0, 5, 10, 15, 20]
    wy = [12, 36]
    plt_title = fr'{mod_label} Uptake: Eugly.'
    reg_mod = 'cmrglc'
else:

    if sys.argv[1] == 'oxy':
        figure = 'figure_3'
        mod = 'oxy_suvr'
        mod_label = r'[$^{15}O$]$O_2$'
        reg_mod = 'om'
    else:
        figure = 'figure_5'
        mod = 'ho_suvr'
        mod_label = r'[$^{15}O$]$H_2O$'
        reg_mod = 'ho'
              
    # Common options for SUVR plots
    cb_label = 'SUVR'
    img_min = 0.5
    img_max = 1.5
    sx =  [0.25, 1.75]
    sy = [-0.1, 0.1]
    sy_ticks = [-.1, -0.05, 0, 0.05, 0.1]
    wy = [0.45, 0.825]
    diff_min = 0
    diff_max = 0.4
    plt_title = fr'{mod_label} Uptake: Eugly.'
    
diff_mid = (diff_max - diff_min) / 2
img_mid = (img_max - img_min) / 2.0 + img_min
x = 48
y = 54
z = 38
dpi = 300.0

#Load in image headers
eugly = nib.load(f'{mod}_basal{suff}.nii.gz')
diff = nib.load(f'{mod}_hypergly_coef{suff}.nii.gz')
mpr = nib.load(f'../common/template_2mm.nii.gz')
mask = nib.load(f'../common/template_2mm_mask.nii.gz')
roi = nib.load(f'../common/MNI152_wmparc_comb_2mm_on_icbm_dwm.nii.gz')
fdr = nib.load(f'{mod}_hypergly_logp{suff}_fdr_05.nii.gz')

# Load in white matter data
reg_data = pd.read_excel('../common/data_values_file.xlsx', 'Freesurfer.Rois')
demo_data = pd.read_excel('../common/data_values_file.xlsx', 'Demographics')
wm_data = reg_data.loc[
    (reg_data.Region == 'Deep White Matter') & (reg_data.Modality == reg_mod)
]
wm_data = pd.merge(
    demo_data, wm_data, how='right', on=['ID', 'Visit.Order', 'Condition'],
)

# Update cmrglc to account for lumped constant if needed
if sys.argv[1] == 'cmrglc_lc':
    lc = -0.0043 * (wm_data['Plasma.Glu'] / 18.0156) + 0.8315
    wm_data['Value'] *= 0.81 / lc

# Average repeat visits and pivot to wide
wm_data = wm_data.groupby(['ID', 'Condition'])['Value'].mean().reset_index()
wm_data = wm_data.pivot(index='ID', columns='Condition', values='Value')

#Load in image data
eugly_data = np.squeeze(eugly.get_fdata())
diff_data = np.squeeze(diff.get_fdata())
fdr_data = np.squeeze(fdr.get_fdata())
mpr_data = np.float32(np.squeeze(mpr.get_fdata()))
mask_data = np.logical_or(np.squeeze(mask.get_fdata()) == 0.0, eugly_data == 0.0)
roi_data = np.float32(np.squeeze(roi.get_fdata()))

# Get average of basal/hyperglycemia
avg_data = eugly_data + diff_data / 2

#Apply brain mask
eugly_data[mask_data] = np.nan
diff_data[mask_data] = np.nan
avg_data[mask_data] = np.nan
fdr_data[mask_data] = np.nan
mpr_data[mask_data] = np.nan
roi_data[mask_data] = np.nan
roi_data[roi_data == 0] = np.nan

#Make colormap for average images
img_map = plt.get_cmap("plasma")
img_map.set_bad("white", alpha=0)
img_map.set_under(img_map(0))
img_map.set_over(img_map(255))

#Make structural colormap
mpr_map = plt.get_cmap("gray")
mpr_map.set_bad("white", alpha=0)

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
ax1.imshow(np.rot90(eugly_data[x, :, :]), cmap=img_map, vmin=img_min, vmax=img_max)
ax2.imshow(np.rot90(mpr_data[:, :, z]), cmap=mpr_map)
ax2.imshow(np.rot90(eugly_data[:, :, z]), cmap=img_map, vmin=img_min, vmax=img_max)
ax3.imshow(np.rot90(mpr_data[:, y, :]), cmap=mpr_map)
ax3.imshow(np.rot90(eugly_data[:, y, :]), cmap=img_map, vmin=img_min, vmax=img_max)

#Make first row
pax1 = ax1.get_position()
pax2 = ax2.get_position()
pax3 = ax3.get_position()
ax1.set_position([0.1, 0.125, pax1.width, pax1.height])
ax2.set_position([0.38, 0.08, pax2.width, pax2.height])
ax3.set_position([0.615, 0.12, pax3.width, pax3.height])

#Add title
ax1.text(0.1, 0.725, 'A)', weight='bold', size=8, transform=fig1.transFigure)
ax2.set_title(plt_title, weight='bold', size=8)

#Add colorbar
cb = mpl.colorbar.ColorbarBase(ax4, cmap=img_map, ticks=[0, 0.5, 1])
cb.ax.set_position([0.875, 0.33, 0.0225, 0.375])
cb.set_ticklabels(['%.1f'%(i) for i in [img_min, img_mid, img_max]])
for tick in cb.ax.yaxis.get_major_ticks():
    tick.label2.set_weight('bold')
    tick.label2.set_size(4.5)
cb.ax.set_xlabel(fr'${cb_label}$', size=5.5, weight='bold', labelpad=5)

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
ax_pos = fig2.add_subplot(sp2[0, 3])
ax_neg = fig2.add_subplot(sp2[0, 3])

# Dilate thresholded image so we can see boundary. Do it for each orthogonal slice
# so we aren't getting anything coming from another slice
diff_thr = np.copy(diff_data)
diff_thr[np.abs(fdr_data) < 1.3] = np.nan
diff_thr_mask = np.float32(~np.isnan(diff_thr))
diff_thr_dil_x = np.float32(binary_dilation(diff_thr_mask[x, :, :])) - diff_thr_mask[x, :, :]
diff_thr_dil_y = np.float32(binary_dilation(diff_thr_mask[:, y, :])) - diff_thr_mask[:, y, :]
diff_thr_dil_z = np.float32(binary_dilation(diff_thr_mask[:, :, z])) - diff_thr_mask[:, :, z]

# Mask out zeros
diff_thr_dil_x[np.logical_or(diff_thr_dil_x == 0, mask_data[x, :, :])] = np.nan
diff_thr_dil_y[np.logical_or(diff_thr_dil_y == 0, mask_data[:, y, :])] = np.nan
diff_thr_dil_z[np.logical_or(diff_thr_dil_z == 0, mask_data[:, :, z])] = np.nan

#Make second row
ax6.imshow(np.rot90(mpr_data[x, :, :]), cmap=mpr_map)
ax6.imshow(
    np.rot90(diff_data[x, :, :]), cmap=pos_map, vmin=diff_min, vmax=diff_max, alpha=0.3
)
ax6.imshow(
    np.rot90(-diff_data[x, :, :]), cmap=neg_map, vmin=diff_min, vmax=diff_max, alpha=0.3
)
ax6.imshow(np.rot90(diff_thr[x, :, :]), cmap=pos_map, vmin=diff_min, vmax=diff_max)
ax6.imshow(np.rot90(-diff_thr[x, :, :]), cmap=neg_map, vmin=diff_min, vmax=diff_max)
ax6.imshow(np.rot90(diff_thr_dil_x), cmap='gray', vmax=1)
ax7.imshow(np.rot90(mpr_data[:, :, z]), cmap=mpr_map)
ax7.imshow(
    np.rot90(diff_data[:, :, z]), cmap=pos_map, vmin=diff_min, vmax=diff_max, alpha=0.3
)
ax7.imshow(
    np.rot90(-diff_data[:, :, z]), cmap=neg_map, vmin=diff_min, vmax=diff_max, alpha=0.3
)
ax7.imshow(np.rot90(diff_thr[:, :, z]), cmap=pos_map, vmin=diff_min, vmax=diff_max)
ax7.imshow(np.rot90(-diff_thr[:, :, z]), cmap=neg_map, vmin=diff_min, vmax=diff_max)
ax7.imshow(np.rot90(diff_thr_dil_z), cmap='gray', vmax=1)
ax8.imshow(np.rot90(mpr_data[:, y, :]), cmap=mpr_map)
ax8.imshow(
    np.rot90(diff_data[:, y, :]), cmap=pos_map, vmin=diff_min, vmax=diff_max, alpha=0.3
)
ax8.imshow(
    np.rot90(-diff_data[:, y, :]), cmap=neg_map, vmin=diff_min, vmax=diff_max, alpha=0.3
)
ax8.imshow(np.rot90(diff_thr[:, y, :]), cmap=pos_map, vmin=diff_min, vmax=diff_max)
ax8.imshow(np.rot90(-diff_thr[:, y, :]), cmap=neg_map, vmin=diff_min, vmax=diff_max)
ax8.imshow(np.rot90(diff_thr_dil_y), cmap='gray', vmax=1)

#Move second row
pax6 = ax6.get_position()
pax7 = ax7.get_position()
pax8 = ax8.get_position()
ax6.set_position([0.1, 0.125, pax6.width, pax6.height])
ax7.set_position([0.38, 0.08, pax7.width, pax7.height])
ax8.set_position([0.615, 0.12, pax8.width, pax8.height])

#Add title
ax6.text(0.1, 0.725, 'B)', weight='bold', size=8, transform=fig2.transFigure)
ax7.set_title('Hygly. - Eugly.', weight='bold', size=8)

#Add positive colorbar
cb2 = mpl.colorbar.ColorbarBase(ax_pos, cmap=pos_map, ticks=[0, 0.5, 1])
cb2.ax.set_position([0.875, 0.57, 0.0225, 0.175])
cb2.set_ticklabels(['% .1f'%(i) for i in [diff_min, diff_mid, diff_max]])
for tick in cb2.ax.yaxis.get_major_ticks():
    tick.label2.set_weight('bold')
    tick.label2.set_size(4.5)

#Add negative colorbar
cb3 = mpl.colorbar.ColorbarBase(ax_neg, cmap=disp_map, ticks=[0, 0.5, 1])
cb3.ax.set_position([0.875, 0.33, 0.0225, 0.175])
cb3.set_ticklabels(['% .1f'%(i) for i in [-diff_max, -diff_mid, -diff_min]])
for tick in cb3.ax.yaxis.get_major_ticks():
    tick.label2.set_weight('bold')
    tick.label2.set_size(4.5)
cb3.ax.set_xlabel(fr'${cb_label}$', size=5.5, weight='bold', labelpad=3)

#Turn off the axis labels
ax6.set_axis_off()
ax7.set_axis_off()
ax8.set_axis_off()

#Save figure
plt.savefig(f'{figure}b.png', dpi=dpi, bbox_inches='tight')
plt.close()

# Make scatter plots to figure
fig3 = plt.figure(figsize=[1500.0 / dpi, 600.0 / dpi])

#Make first grid
sp3 = gridspec.GridSpec(ncols=2, nrows=1, figure=fig3, wspace=0.325)

#Add subplots to grid
ax10 = fig3.add_subplot(sp3[0, 0])
ax11 = fig3.add_subplot(sp3[0, 1])

# Make white matter spagetti plot
ax10.grid(linewidth=0.25, c='gray', which='major')
for index, row in wm_data.iterrows():
    ax10.plot(row.index, row.values, marker='o', c='#b6006b', ms=1.25, lw=1)
ax10.set_xlabel('Condition', fontweight='bold', fontsize=5)
ax10.set_xlim([-0.3, 1.3])
ax10.set_xticks(['basal', 'hypergly'], ['Eugly.', 'Hygly.'])
ax10.set_ylabel(fr'{mod_label} (${cb_label}$)', fontweight='bold', fontsize=5)
ax10.set_ylim(wy)
ax10.set_title('Deep White Matter', fontweight='bold', fontsize=5.5, y=0.985)
ax10.tick_params(axis='both', labelsize=4.5, width=0.5, colors='gray')
for axis in ['top','bottom','left','right']:
    ax10.spines[axis].set_linewidth(0.5)
    ax10.spines[axis].set_color('gray')
ax10.text(-0.225, 1.1, 'C)', weight='bold', size=8, transform=ax10.transAxes)

# Add white matter image figure
roi_ax = ax10.inset_axes([.05, .575, .5, .5])
roi_ax.matshow(np.rot90(mpr_data[:, y, :]), cmap=mpr_map)
roi_ax.matshow(np.rot90(roi_data[:, y, :]), cmap=roi_map)
roi_ax.set_axis_off()

# Compute density of baseline vs. change 
dist = np.stack((avg_data.flatten(), diff_data.flatten()), axis=1)
dist = dist[np.isfinite(dist[:, 0]), :]
pdf = fastkde.pdf(
    dist[:, 0], dist[:, 1], var_names = ['avg', 'diff'], num_points=2**8 + 1
)

# Interpolate pdf so that each voxel gets a value
pdf_vals = pdf.to_numpy().T
pdf_i = interpn(
    (pdf.coords['avg'].data, pdf.coords['diff'].data),
    pdf_vals, dist,
    bounds_error=False,
    method='linear',
    fill_value=None
)
pdf_i /= np.max(pdf_i)

# Correlate baseline and change
ax11.grid(linewidth=0.25, c='gray')
ax11.scatter(dist[:, 0], dist[:, 1], c=pdf_i, alpha=pdf_i, s=0.5, cmap='plasma', zorder=2)
ax11.set_xlabel(rf'Avg. {mod_label} (${cb_label}$)', fontweight='bold', fontsize=5)
ax11.set_ylabel(
    rf'$\Delta$ {mod_label} (${cb_label}$)', fontweight='bold', fontsize=5, labelpad=0
)
ax11.set_xlim(sx)
ax11.set_ylim(sy)
ax11.set_yticks(sy_ticks)
ax11.set_title(rf'{mod_label} Spatial Correlation', fontweight='bold', fontsize=5.5)
ax11.tick_params(axis='both', labelsize=4.5, labelcolor='gray')
ax11.tick_params(axis='both', labelsize=4.5, width=0.5, colors='gray')
for axis in ['top','bottom','left','right']:
    ax11.spines[axis].set_linewidth(0.5)
    ax11.spines[axis].set_color('gray')
ax11.text(-0.225, 1.1, 'D)', weight='bold', size=8, transform=ax11.transAxes)
cor = np.corrcoef(dist[:, 0], dist[:, 1])[0, 1]
ax11.text(0.7, 0.8, f'r = {cor:.2f}', weight='bold', size=6, transform=ax11.transAxes)

#Save figure
plt.savefig(f'{figure}cd.png', dpi=dpi, bbox_inches='tight')
plt.close('all')

#Read in each image
hyper_img = Image.open(f'{figure}a.png')
eugly_img = Image.open(f'{figure}b.png')
plot_img = Image.open(f'{figure}cd.png')

#Crop voxelwise images
coords = (0, 0, hyper_img.size[0], hyper_img.size[1] * 0.9)
hygly_crop = hyper_img.crop(coords)
eugly_crop = eugly_img.crop(coords)

#Resample the scatter plot figure
plot_resamp = plot_img.resize(
    (int(float(plot_img.size[0]*0.95)), int(float(plot_img.size[1]*0.95))), Image.LANCZOS
)

#Make new image
output = Image.new("RGBA",
                   (hygly_crop.size[0], int(hygly_crop.size[1] * 3.4)),
                   color=(255, 255, 255))

#Paste images to output
output.paste(hygly_crop, (0, 0))
output.paste(eugly_crop, (0, eugly_crop.size[1]))
output.paste(plot_resamp, (10, int(eugly_crop.size[1] * 2.025)))

#Save image
output.save(f'{figure}.tiff')
for f in [f'{figure}a.png', f'{figure}b.png', f'{figure}cd.png']:
    os.remove(f)