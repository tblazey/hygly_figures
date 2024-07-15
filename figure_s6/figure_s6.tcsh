#!/bin/tcsh -f
echo $cwd


python3 ../common/three_panel_scatter.py cbf
Rscript ../common/wm_summary.R cbf
