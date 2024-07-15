import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
from PIL import Image

#Change default math font
mpl.rcParams['mathtext.default'] = 'regular'

#Define function for computing bound
def calc_bound(img_data, max_p):

    #Get minimum and maximum at range
    pos_bound = np.nanpercentile(img_data[img_data > 0], max_p)
    neg_bound = np.nanpercentile(img_data[img_data < 0], max_p) * -1.0
    bound = np.max((pos_bound, neg_bound))
    
    return bound

#Set coordinates for ploting
x = 48
y = 54
z = 38

#Load in image headers
quant =  nib.load('../figure_2/cmrglc_hypergly_coef.nii.gz')
suvr_wb = nib.load('fdg_hypergly_coef_wb_no_smooth.nii.gz')
suvr_po = nib.load('fdg_hypergly_coef_no_smooth.nii.gz')
mask = nib.load('../common/template_2mm_mask.nii.gz')
roi = nib.load('MNI152_wmparc_comb_2mm_on_icbm_po.nii.gz')
mpr = nib.load('../common/template_2mm.nii.gz')

#Load in image data
quant_data = np.squeeze(quant.get_fdata())
suvr_wb_data = np.squeeze(suvr_wb.get_fdata())
suvr_po_data = np.squeeze(suvr_po.get_fdata())
mask_data = np.squeeze(mask.get_fdata()) == 0.0
roi_data = np.float32(np.squeeze(roi.get_fdata()))
mpr_data = np.squeeze(mpr.get_fdata())

#Apply mask
quant_data[mask_data] = np.nan
suvr_wb_data[mask_data] = np.nan
suvr_po_data[mask_data] = np.nan
roi_data[mask_data] = np.nan
roi_data[roi_data == 0.0] = np.nan
mpr_data[mask_data] = np.nan

#Get bound for images
quant_bound = np.nanpercentile(quant_data, 98)
wb_bound = np.nanpercentile(suvr_wb_data, 98)
po_bound = np.nanpercentile(suvr_po_data, 98)

#Make figure object
dpi = 300.0
fig1 = plt.figure(figsize=[2000 / dpi, 800.0 / dpi])

#Make first grid
sp1 = gridspec.GridSpec(ncols=4, nrows=1, figure=fig1,
                        width_ratios=[1, 1, 1, 0.1])

#Add subplots to grid
ax1 = fig1.add_subplot(sp1[0, 0])
ax2 = fig1.add_subplot(sp1[0, 1])
ax3 = fig1.add_subplot(sp1[0, 2])
ax4 = fig1.add_subplot(sp1[0, 3])

#Show images
ax1.imshow(np.rot90(quant_data[:, :, z]), cmap='RdBu_r', vmin=-quant_bound, vmax=quant_bound)
ax2.imshow(np.rot90(suvr_wb_data[:, :, z]), cmap='RdBu_r', vmin=-wb_bound, vmax=wb_bound)
ax3.imshow(np.rot90(suvr_po_data[:, :, z]), cmap='RdBu_r', vmin=-po_bound, vmax=po_bound)

#Move images
pax1 = ax1.get_position()
pax2 = ax2.get_position()
pax3 = ax3.get_position()
ax1.set_position([0.025, 0.10, pax1.width, pax1.height])
ax2.set_position([0.275, 0.10, pax2.width, pax2.height])
ax3.set_position([0.525, 0.10, pax3.width, pax3.height])

#Add image labels
ax1.text(0.009, 0.99, 'C)', weight='bold', size=8, transform=ax1.transAxes)
ax1.text(0.30, 0.99, 'CMRglc', weight='bold', size=8, transform=ax1.transAxes)
ax1.text(0.009, 0.99, 'D)', weight='bold', size=8, transform=ax2.transAxes)
ax1.text(0.30, 0.99, r'$SUVR_{wb}$', weight='bold', size=8, transform=ax2.transAxes)
ax1.text(0.009, 0.99, 'E)', weight='bold', size=8, transform=ax3.transAxes)
ax1.text(0.30, 0.99, r'$SUVR_{ref}$', weight='bold', size=8, transform=ax3.transAxes)

#Add colorbar
cb1 = mpl.colorbar.ColorbarBase(ax4, cmap='RdBu_r', ticks=[0, 0.5, 1])
cb1.ax.set_position([0.775, 0.325, 0.0175, 0.275])
cb1.set_ticklabels(['%s'%(i) for i in ['-98 %', 0, '+98%']])
for tick in cb1.ax.yaxis.get_major_ticks():
    tick.label2.set_size(5.0)

#Turn off the axis labels
ax1.set_axis_off()
ax2.set_axis_off()
ax3.set_axis_off()

#Save figure
plt.savefig('figure_s2cde.png', dpi=dpi, bbox_inches='tight')
plt.close('all')

# Make second figure
fig2 = plt.figure(figsize=[2000 / dpi, 1000.0 / dpi])

#Make first grid
sp2 = gridspec.GridSpec(ncols=2, nrows=1, figure=fig2, width_ratios=[1, 0.7])

#Add subplots to grid
ax5 = fig2.add_subplot(sp2[0, 0])
ax6 = fig2.add_subplot(sp2[0, 1])

# Load in images with roi values
roi_means = nib.load('cmrglc_wmparc_hypergly_coef.nii.gz').get_fdata().flatten()
roi_ses = nib.load('cmrglc_wmparc_hypergly_se.nii.gz').get_fdata().flatten()
roi_names = pd.read_excel('../common/data_values_file.xlsx', "Allen.Gene", usecols=[5])

# Add whole-brain data
wb_df = pd.read_csv('../table_s3/table_s3.csv', index_col=0)
roi_means = np.append(roi_means, wb_df.Diff.CMRglc)
roi_cis = np.append(roi_ses * 1.96, wb_df['Diff.CI'].CMRglc)
roi_names = np.append(roi_names, 'Whole Brain')

# Create color code
cols = np.repeat('#006bb6', roi_names.shape[0])
cols[-1] = '#000000'
cols[np.argsort(np.abs(roi_means))[0:3]] = '#6bb600'

# Make plot
roi_sort = np.argsort(roi_means)
ax5.grid()
ax5.grid(linewidth=0.25, c='gray')
ax5.scatter(
    roi_means[roi_sort], roi_names[roi_sort], s=4, zorder=3, c=cols[roi_sort]
)
ax5.errorbar(
    roi_means[roi_sort],
    roi_names[roi_sort],
    xerr=roi_cis[roi_sort],
    zorder=2,
    lw=0.75,
    c='black',
    ls='none'
)
ax5.set_xlabel(
    r'$\Delta$CMRglc ($\mu$Mol/hg/min)', fontweight='bold', fontsize=5
)
ax5.set_title('Regional CMRglc', fontweight='bold', fontsize=7)
ax5.tick_params(
    axis='both',
    labelsize=4.5,
    labelcolor='gray',
    width=0.25,
    colors='gray')
ax5.tick_params(axis='y', labelsize=4)
for axis in ['top','bottom','left','right']:
    ax5.spines[axis].set_linewidth(0.5)
    ax5.spines[axis].set_color('gray')
ax5.text(-0.3, 1.05, 'A)', weight='bold', size=8, transform=ax5.transAxes)

#Make colormap for roi
roi_map = mpl.colors.ListedColormap([(0.42, 0.714, 0.0)])

#Make structural colormap
mpr_map = plt.get_cmap("gray")
mpr_map.set_bad("white", alpha=0)

# Add coronal image
cor_ax = ax5.inset_axes([.675, .3, .3, .3])
cor_ax.matshow(np.rot90(mpr_data[:, 59, :]), cmap=mpr_map)
cor_ax.matshow(np.rot90(roi_data[:, 59, :]), cmap=roi_map)
cor_ax.set_axis_off()

# Add sagittal image
sag_ax = ax5.inset_axes([.645, .075, .35, .35])
sag_ax.matshow(np.rot90(mpr_data[66, :, :]), cmap=mpr_map)
sag_ax.matshow(np.rot90(roi_data[66, :, :]), cmap=roi_map)
sag_ax.set_axis_off()

# Load in reference region data
ref_df = pd.read_excel('../common/data_values_file.xlsx', "Ref.CMRglc")

# Make scatter plot for reference CMRglc
ax6.grid(linewidth=0.25, c='gray')
for _, group in ref_df.groupby(['ID']):
    sub_df = group.groupby(['ID', 'Condition']).mean().reset_index()
    ax6.plot(
        sub_df.Condition, sub_df.Value, c='#6bb600', marker='o', ms=2, lw=1, zorder=2
    )

# Format plot
ax6.set_ylabel(
    r'CMRglc ($\mu$Mol/hg/min)', fontweight='bold', fontsize=5
)
ax6.set_xlabel('Condition', fontweight='bold', fontsize=5)
ax6.set_xticks(['basal', 'hypergly'], labels=['Eugly.', 'Hygly.'])
ax6.set_xlim([-0.3, 1.3])
ax6.invert_xaxis()
ax6.set_title('Reference Region CMRglc', fontweight='bold', fontsize=7)
ax6.tick_params(
    axis='both',
    labelsize=4.5,
    labelcolor='gray',
    width=0.25,
    colors='gray')
for axis in ['top','bottom','left','right']:
    ax6.spines[axis].set_linewidth(0.5)
    ax6.spines[axis].set_color('gray')
ax6.text(-0.2, 1.05, 'B)', weight='bold', size=8, transform=ax6.transAxes)

# Then save it
plt.subplots_adjust(wspace=0.3)
plt.savefig('figure_s2ab.png', dpi=dpi, bbox_inches='tight')
plt.close('all')

#Load in roi image
top_img = Image.open('figure_s2ab.png')
bot_img = Image.open('figure_s2cde.png')
                        
#Make container for pasting images into
output = Image.new(
    "RGBA",
    (top_img.width, top_img.height + bot_img.height),
    color=(255, 255, 255)
)

#Paste images
output.paste(top_img, (0, 0))
output.paste(bot_img, (100, top_img.height + 50))

#Save output image
output.save('figure_s2.tiff')






