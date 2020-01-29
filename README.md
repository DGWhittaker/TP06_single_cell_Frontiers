# TP06 single cell code

This repository contains C++ code for [TP06 single cell](https://doi.org/10.1152/ajpheart.00109.2006) models and Python plotting scipts used in our [recent Frontiers paper](https://doi.org/10.3389/fphys.2019.00308)

* Compile with single cell model with `g++ TNNP.cpp -o ttcell`. Short and basal action potential variants are described within code.
* Run `python plot_APs.py` or `python plot_APs_short.py` to output figures.
* [APs](https://github.com/DGWhittaker/TP06_single_cell_Frontiers/tree/master/APs) contains single cell model output text files required for figures.

# Acknowledging this work

If you publish any work based on the contents of this repository please cite:

Perez Alday, E. A., Whittaker, D. G., Benson, A. P., Colman, M. A.
(2019).
[Effects of Heart Rate and Ventricular Wall Thickness on Non-invasive Mapping: An _in silico_ Study](https://doi.org/10.3389/fphys.2019.00308).
_Frontiers in Physiology_, 10, 308.