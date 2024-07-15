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

if sys.argv[1] == 'cmrglc':
    figure = 'figure_2'
    suff = ''
elif sys.argv[1] == 'cmrglc_lc':
    figure = 'figure_s1'
    suff = '_lc'
else:
    sys.exit()

#Load in image headers
mpr = nib.load('../common/template_2mm.nii.gz')
mask = nib.load(f'../common/template_2mm_mask.nii.gz')
roi = nib.load(f'../common/MNI152_wmparc_comb_2mm_on_icbm_dwm.nii.gz')
eugly = nib.load(f'cmrglc_basal{suff}.nii.gz')
diff = nib.load(f'cmrglc_hypergly_coef{suff}.nii.gz')
fdr = nib.load(f'cmrglc_hypergly_logp{suff}_fdr_05.nii.gz')

# Get regional data and demographics
reg_data = pd.read_excel('../common/data_values_file.xlsx', 'Freesurfer.Rois')
demo_data = pd.read_excel('../common/data_values_file.xlsx', 'Demographics')
wm_data = reg_data.loc[
    (reg_data.Region == 'Deep White Matter') & (reg_data.Modality == 'cmrglc')
]
wm_data = pd.merge(
    demo_data, wm_data, how='right', on=['ID', 'Visit.Order', 'Condition'],
)

# Update cmrglc to account for lumped constant if needed
if sys.argv[1] == 'cmrglc_lc':
    lc = -0.0043 * (wm_data['Plasma.Glu'] / 18.0156) + 0.8315
    wm_data['Value'] *= 0.81 / lc

# Average if a participant has multiple values 
wm_data = wm_data.groupby(['ID', 'Condition']).mean(['Condition']).reset_index()

#Set coordinates for ploting
x = 48
y = 54
z = 38

#Set range for image
eugly_min = 10
eugly_max = 50
eugly_mid = (eugly_max - eugly_min) / 2.0 + eugly_min
diff_min = 2
diff_max = 15
diff_mid = (diff_max - diff_min) / 2.0 + diff_min

#Load in image data
eugly_data = np.squeeze(eugly.get_fdata())
diff_data = np.squeeze(diff.get_fdata())
fdr_data = np.squeeze(fdr.get_fdata())
mpr_data = np.float32(np.squeeze(mpr.get_fdata()))
mask_data = np.logical_or(np.squeeze(mask.get_fdata()) == 0.0, eugly_data == 0.0)
roi_data = np.float32(np.squeeze(roi.get_fdata()))

#Mask outside of fdr threshold
diff_thr = np.copy(diff_data)
diff_thr[np.abs(fdr_data) < 1.3] = np.nan

#Apply brain mask
eugly_data[mask_data] = np.nan
diff_data[mask_data] = np.nan
diff_thr[mask_data] = np.nan
fdr_data[mask_data] = np.nan
mpr_data[mask_data] = np.nan
roi_data[mask_data] = np.nan

roi_data[roi_data == 0.0] = np.nan

#Make colormap for average images
eugly_map = plt.get_cmap("plasma")
eugly_map.set_bad("white", alpha=0)
eugly_map.set_under(eugly_map(0))
eugly_map.set_over(eugly_map(255))

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
dpi = 300.0
fig1 = plt.figure(figsize=[1500.0 / dpi, 600.0 / dpi])

#Make first grid
sp1 = gridspec.GridSpec(ncols=5, nrows=1, figure=fig1,
                        width_ratios=[218.0/182.0, 1, 1, 0.05, 0.05], wspace=0.05)

#Add subplots to grid
ax1 = fig1.add_subplot(sp1[0, 0])
ax2 = fig1.add_subplot(sp1[0, 1])
ax3 = fig1.add_subplot(sp1[0, 2])
ax4 = fig1.add_subplot(sp1[0, 3])
ax5 = fig1.add_subplot(sp1[0, 4])

#Make first row
ax1.imshow(np.rot90(mpr_data[x, :, :]), cmap=mpr_map)
ax1.imshow(np.rot90(diff_thr[x, :, :]), cmap=pos_map, vmin=diff_min, vmax=diff_max)
ax1.imshow(np.rot90(-diff_thr[x, :, :]), cmap=neg_map, vmin=diff_min, vmax=diff_max)
ax2.imshow(np.rot90(mpr_data[:, :, z]), cmap=mpr_map)
ax2.imshow(np.rot90(diff_thr[:, :, z]), cmap=pos_map, vmin=diff_min, vmax=diff_max)
ax2.imshow(np.rot90(-diff_thr[:, :, z]), cmap=neg_map, vmin=diff_min, vmax=diff_max)
ax3.imshow(np.rot90(mpr_data[:, y, :]), cmap=mpr_map)
ax3.imshow(np.rot90(diff_thr[:, y, :]), cmap=pos_map, vmin=diff_min, vmax=diff_max)
ax3.imshow(np.rot90(diff_thr[:, y, :]), cmap=pos_map, vmin=diff_min, vmax=diff_max)
ax3.imshow(np.rot90(-diff_thr[:, y, :]), cmap=neg_map, vmin=diff_min, vmax=diff_max)

#Make first row
pax1 = ax1.get_position()
pax2 = ax2.get_position()
pax3 = ax3.get_position()
ax1.set_position([0.1, 0.125, pax1.width, pax1.height])
ax2.set_position([0.38, 0.08, pax2.width, pax2.height])
ax3.set_position([0.615, 0.12, pax3.width, pax3.height])

#Add title
ax1.text(0.1, 0.725, 'B)', weight='bold', size=8, transform=fig1.transFigure)
ax2.set_title('Hyper. - Eugly.', weight='bold', size=8)

#Add positive colorbar
cb3 = mpl.colorbar.ColorbarBase(ax4, cmap=pos_map, ticks=[0, 0.5, 1])
cb3.ax.set_position([0.875, 0.57, 0.0225, 0.175])
cb3.set_ticklabels(['% .1f'%(i) for i in [diff_min, diff_mid, diff_max]])
for tick in cb3.ax.yaxis.get_major_ticks():
    tick.label2.set_weight('bold')
    tick.label2.set_size(4.5)

#Add negative colorbar
cb4 = mpl.colorbar.ColorbarBase(ax5, cmap=disp_map, ticks=[0, 0.5, 1])
cb4.ax.set_position([0.875, 0.33, 0.0225, 0.175])
cb4.set_ticklabels(['% .1f'%(i) for i in [-diff_max, -diff_mid, -diff_min]])
for tick in cb4.ax.yaxis.get_major_ticks():
    tick.label2.set_weight('bold')
    tick.label2.set_size(4.5)
cb_label = '\\frac{\mu Mol}{hg \cdot min}'
cb4.ax.set_xlabel(r'$%s$'%(cb_label), size=5.5, weight='bold', labelpad=3)

#Turn off the axis labels
ax1.set_axis_off()
ax2.set_axis_off()
ax3.set_axis_off()

#Save figure
plt.savefig(f'{figure}b.png', dpi=dpi, bbox_inches='tight')
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
ax6.imshow(np.rot90(eugly_data[x, :, :]), cmap=eugly_map, vmin=eugly_min, vmax=eugly_max)
ax7.imshow(np.rot90(mpr_data[:, :, z]), cmap=mpr_map)
ax7.imshow(np.rot90(eugly_data[:, :, z]), cmap=eugly_map, vmin=eugly_min, vmax=eugly_max)
ax8.imshow(np.rot90(mpr_data[:, y, :]), cmap=mpr_map)
ax8.imshow(np.rot90(eugly_data[:, y, :]), cmap=eugly_map, vmin=eugly_min, vmax=eugly_max)

#Make second row
pax6 = ax6.get_position()
pax7 = ax7.get_position()
pax8 = ax8.get_position()
ax6.set_position([0.1, 0.125, pax6.width, pax6.height])
ax7.set_position([0.38, 0.08, pax7.width, pax7.height])
ax8.set_position([0.615, 0.12, pax8.width, pax8.height])

#Add title
ax6.text(0.1, 0.725, 'A)', weight='bold', size=8, transform=fig2.transFigure)
ax7.set_title('CMRglc: Eugly.', weight='bold', size=8)

#Add colorbar
cb2 = mpl.colorbar.ColorbarBase(ax9, cmap=eugly_map, ticks=[0, 0.5, 1])
cb2.ax.set_position([0.875, 0.33, 0.0225, 0.375])
cb2.set_ticklabels(['%.1f'%(i) for i in [eugly_min, eugly_mid, eugly_max]])
for tick in cb2.ax.yaxis.get_major_ticks():
    tick.label2.set_weight('bold')
    tick.label2.set_size(4.5)
cb2.ax.set_xlabel(r'$%s$'%(cb_label), size=5.5, weight='bold', labelpad=5)

#Turn off the axis labels
ax6.set_axis_off()
ax7.set_axis_off()
ax8.set_axis_off()

#Save figure
plt.savefig(f'{figure}a.png', dpi=dpi, bbox_inches='tight')
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
ax10.set_xticklabels(['Hyper.', 'Eugly.'])
ax10.set_xlim([-0.3, 1.3])
ax10.invert_xaxis()
ax10.set_ylabel(r'CMRglc ($\mu$Mol/hg/min)', fontweight='bold', fontsize=5)
ax10.set_yticks([15, 20, 25, 30, 35])
ax10.set_title('Deep White Matter', fontweight='bold', fontsize=5.5,  y=0.99)
ax10.tick_params(axis='both', labelsize=4.5, width=0.5, colors='gray')
for axis in ['top','bottom','left','right']:
    ax10.spines[axis].set_linewidth(0.5)
    ax10.spines[axis].set_color('gray')
ax10.text(-0.225, 1.1, 'C)', weight='bold', size=8, transform=ax10.transAxes)

# Add white matter image figure
roi_ax = ax10.inset_axes([.1, .55, .5, .5])
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
ax11.scatter(dist[:, 0], dist[:, 1], c=pdf_i, s=0.5, cmap='plasma', alpha=pdf_i, zorder=2)
ax11.set_xlabel(r'Eugly. CMRglc ($\mu$Mol/hg/min)', fontweight='bold', fontsize=5)
ax11.set_ylabel(
    r'$\Delta$CMRglc ($\mu$Mol/hg/min)', fontweight='bold', fontsize=5, labelpad=0
)
ax11.set_title('CMRglc Spatial Correlation', fontweight='bold', fontsize=5.5)
ax11.tick_params(axis='both', labelsize=4.5, labelcolor='gray')
ax11.tick_params(axis='both', labelsize=4.5, width=0.5, colors='gray')
for axis in ['top','bottom','left','right']:
    ax11.spines[axis].set_linewidth(0.5)
    ax11.spines[axis].set_color('gray')
ax11.text(-0.225, 1.1, 'D)', weight='bold', size=8, transform=ax11.transAxes)
cor = np.corrcoef(dist[:, 0], dist[:, 1])[0, 1]
ax11.text(0.6, 0.8, f'r = {cor:.2}', weight='bold', size=6, transform=ax11.transAxes)

#Save figure
plt.savefig(f'{figure}cd.png', dpi=dpi, bbox_inches='tight')
plt.close('all')

#Read in each image
eug_img = Image.open(f'{figure}a.png')
diff_img = Image.open(f'{figure}b.png')
plot_img = Image.open(f'{figure}cd.png')

#Crop voxelwise images
coords = (0, 0, diff_img.size[0], diff_img.size[1] * 0.9)
diff_crop = diff_img.crop(coords)
eug_crop = eug_img.crop(coords)

#Resample the scatter plot figure
plot_resamp = plot_img.resize(
    (int(float(plot_img.size[0]*0.95)), int(float(plot_img.size[1]*0.95))), Image.LANCZOS
)

#Make new image
output = Image.new("RGBA",
                   (diff_crop.size[0], int(diff_crop.size[1] * 3.4)),
                   color=(255, 255, 255))

#Paste images to output
output.paste(eug_crop, (0, 0))
output.paste(diff_crop, (0, eug_crop.size[1]))
output.paste(plot_resamp, (10, int(eug_crop.size[1] * 2.025)))

#Save image
output.save(f'{figure}.tiff')