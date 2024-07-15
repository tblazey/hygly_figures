#!/bin/tcsh -f
echo $cwd

python3 figure_s2.py
Rscript fit_ref.R
