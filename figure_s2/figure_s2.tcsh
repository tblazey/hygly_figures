#!/bin/tcsh -f
echo $cwd

Rscript fit_regions.R
python3 figure_s2.py
Rscript fit_ref.R
