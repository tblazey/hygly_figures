#!/bin/tcsh -f
echo $cwd


python3 ../common/three_panel_cmrglc.py cmrglc
Rscript ../common/wm_summary.R cmrglc
