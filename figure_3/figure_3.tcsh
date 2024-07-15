#!/bin/tcsh -f
echo $cwd


python3 ../common/three_panel_scatter.py oxy
Rscript ../common/wm_summary.R oxy
