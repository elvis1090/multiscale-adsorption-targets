# multiscale-adsorption-targets

## Material property targets for emerging nanomaterials to enable point-of-use and point-of-entry water treatment systems
### Elvis A. Eugene, William A. Phillip, Alexander W. Dowling
#### Department of Chemical and Biomolecular Engineering
#### University of Notre Dame, Notre Dame, IN 46556


## Primary analysis of results in the paper:

The primary code for this analysis is written using MATLAB. 
The execution of the scripts is initiated by running run_adsorption_analysis.m
If only part of the analysis needs to be run, the user may do so by modifying the 
variable 'cases' in run_adsorption_analysis.m 
Extensive comments are provided in the scripts for Case definitions and their selection. 
Comments are also provided in all the functions to enable customization.


## Capacity upper bound curves

The capcity upper bound curves were calculated using an optimization problem modeled using 
the Pyomo package in Python. The scripts for capcity upper bound calculation are pyo_q_ub_curve_init.py 
and pyo_q_ub_curve_opt.py. For accurate results, the user is recommended to run the scripts
using Ipopt solver with ma27 linear solver from the HSL Mathematical Software Library.


## Plots
All plot for the paper were made using Matplotlib in Python. Once the data for the plots are 
generated using the MATLAB and Python scripts as described above, the Jupyter notebook 
adsorption-paper-plots.ipynb may be used to readily plot and view the results of the analysis.