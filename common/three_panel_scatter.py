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
    
if sys.argv[1] == 'cbf':
    cb_label = r'$\frac{mL}{hg \cdot min}$'
    mod = 'cbf'
    mod_label = 'CBF'
    figure = 'figure_s6'
    img_min = 10
    img_max = 100 
    sx =  [0, 110]
    sy = [-25, 25]
    sy_ticks = [-20, -10, 0, 10, 20]
    wy = [5, 50]
    plt_title = rf'CBF: Hyper.'
    reg_mod = 'cbf'

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
    plt_title = fr'{mod_label} Uptake: Hyper'
    

img_mid = (img_max - img_min) / 2.0 + img_min
x = 48
y = 54
z = 38
dpi = 300.0

#Load in image headers
eugly = nib.load(f'{mod}_basal.nii.gz')
hygly = nib.load(f'{mod}_hypergly.nii.gz')
diff = nib.load(f'{mod}_hypergly_coef.nii.gz')
mpr = nib.load(f'../common/template_2mm.nii.gz')
mask = nib.load(f'../common/template_2mm_mask.nii.gz')
roi = nib.load(f'../common/MNI152_wmparc_comb_2mm_on_icbm_dwm.nii.gz')

# Load in white matter data
reg_data = pd.read_excel('../common/data_values_file.xlsx', 'Freesurfer.Rois')
demo_data = pd.read_excel('../common/data_values_file.xlsx', 'Demographics')
wm_data = reg_data.loc[
    (reg_data.Region == 'Deep White Matter') & (reg_data.Modality == reg_mod)
]
wm_data = pd.merge(
    demo_data, wm_data, how='right', on=['ID', 'Visit.Order', 'Condition'],
)

#Load in image data
eugly_data = np.squeeze(eugly.get_fdata())
hygly_data = np.squeeze(hygly.get_fdata())
diff_data = np.squeeze(diff.get_fdata())
mpr_data = np.float32(np.squeeze(mpr.get_fdata()))
mask_data = np.logical_or(np.squeeze(mask.get_fdata()) == 0.0, eugly_data == 0.0)
roi_data = np.float32(np.squeeze(roi.get_fdata()))

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
ax1.text(0.1, 0.725, 'A)', weight='bold', size=8, transform=fig1.transFigure)
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
ax6.text(0.1, 0.725, 'B)', weight='bold', size=8, transform=fig2.transFigure)
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

# Make scatter plots to figure
fig3 = plt.figure(figsize=[1500.0 / dpi, 600.0 / dpi])

#Make first grid
sp3 = gridspec.GridSpec(ncols=2, nrows=1, figure=fig3, wspace=0.325)

#Add subplots to grid
ax10 = fig3.add_subplot(sp3[0, 0])
ax11 = fig3.add_subplot(sp3[0, 1])

# Make white matter spagetti plot
ax10.grid(linewidth=0.25, c='gray', which='major')
u_sub = np.unique(wm_data['ID'])
for sub in u_sub:
    sub_mask = wm_data['ID'] == sub
    
    ax10.plot(
        wm_data['Condition'][sub_mask],
        wm_data['Value'][sub_mask],
        marker='o',
        c='#b6006b',
        ms=1.25,
        lw=1
    )
ax10.set_xlabel('Condition', fontweight='bold', fontsize=5)
ax10.set_xticklabels(['Eugly.', 'Hygly.'])
ax10.set_xlim([-0.3, 1.3])
ax10.set_ylabel(fr'{mod_label} ({cb_label})', fontweight='bold', fontsize=5)
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
dist = np.stack((eugly_data.flatten(), diff_data.flatten()), axis=1)
dist = dist[np.isfinite(dist[:, 0]), :]
pdf = fastkde.pdf(
    dist[:, 0], dist[:, 1], var_names = ['eugly', 'diff'], num_points=2**8 + 1
)

# Interpolate pdf so that each voxel gets a value
pdf_vals = pdf.to_numpy().T
pdf_i = interpn(
    (pdf.coords['eugly'].data, pdf.coords['diff'].data),
    pdf_vals, dist,
    bounds_error=False,
    method='linear',
    fill_value=None
)
pdf_i /= np.max(pdf_i)

# Correlate baseline and change
ax11.grid(linewidth=0.25, c='gray')
ax11.scatter(dist[:, 0], dist[:, 1], c=pdf_i, alpha=pdf_i, s=0.5, cmap='plasma', zorder=2)
ax11.set_xlabel(rf'Eugly. {mod_label} ({cb_label})', fontweight='bold', fontsize=5)
ax11.set_ylabel(
    rf'$\Delta$ {mod_label} ({cb_label})', fontweight='bold', fontsize=5, labelpad=0
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