GREAT3 metric simulation and evaluation scripts
===============================================

The scripts in this directory were used to simulate the behaviour of metrics for
the GREAT3 challenge, and to evaluate them for the GREAT3 Leaderboard web page.

## Metric evaluation

The file `evaluate.py` contains the code used to calculate metrics for GREAT3.
It requires the full set of GREAT3 truth data, available from
http://great3.projects.phys.ucl.ac.uk/leaderboard/data/public/great3_metric_evaluation_truth.tar

If you wish to calculate the metrics yourselves please see the `q_constant()`
and `q_variable()` functions of that module, and their docstrings.

## Metric simulations

The following scripts were used to tabulate the bias response of the final
versions of the GREAT3 metrics, following the revision of February 2014:

    tabulate_const_shear_metric_rev1.py
    tabulate_correlated_const_shear_metric_rev1.py
    tabulate_correlated_variable_shear_metric_rev1.py
    tabulate_variable_shear_metric_rev1.py

These use the GREAT3 truth data to run simulations of metric performance for
submissions of varying (linear, m & c-style) bias.  These, and almost all the
metric simulation scripts, make use of utility functions in the module
`g3metrics.py`, also in this directory.

The scripts above with "correlated" in the filename make use of outputs from the
following scripts:

    calculate_const_i3_rg_correlation.py
    calculate_variable_i3_rg_correlation.py

These are used to calculate the product moment correlation coefficient between
submissions from different methods (im3shape and the regauss example script).
For the calculation of variable shear metric simulations with correlations,
results from these calculations are used in

    calculate_variable_covariance_matrix.py
    calculate_variable_cholesky.py

There are also a number of older scripts used to investigate a number of
potential metrics during the preparation for GREAT3, in approximate
chronological order:

    calculate_plot_const_shear_metrics.py
    plot_const_shear_metrics_2D.py
    calculate_QG10_var_shear_normalization.py
    calculate_QG3S_var_shear_normalization.py
    plot_PS_normalization_results.py
    test_evaluate.py
    plot_test_evaluate_metrics.py

Please see the individual module docstrings for more information about these
older simulations.
