# MDDC 1.1.0

In this version, we have updated the output for the the `mddc_mc()` method.

We now include additional output element: 

- `fisher_pval`: p-values for each cell in the step 2 of the algorithm, 
calculated using the Monte Carlo method for cells with count greater than five 
and p-values obtained from the Fisher’s exact test for cells with count less 
than or equal to five in the contingency 
table.

Another output element `mc_pval` is now updated and and now outputs: p-values 
for each cell in the step 2 of the algorithm, calculated using the Monte Carlo 
method.

## New features

# MDDC 1.0.0

This is the first release of MDDC. 
