#!/bin/tcsh -f
echo $cwd


python3 ../common/three_panel_scatter.py cmrglc
Rscript ../common/wm_summary.R cmrglc
