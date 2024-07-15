#!/bin/tcsh -f
echo $cwd


python3 ../common/three_panel_scatter.py ho
Rscript ../common/wm_summary.R ho
