# Interim Augmented Inverse Probability Weighted Estimator

This repository holds files related to sequential analysis using the AIPW estimator in SMARTs. 

The folder 'smartim' can be built as an R package. This library can then be used for simulations. It is depricated. 
The functions in it do not account for interim feasible sets functions. 

The folder 'value_estimator' contains the scripts required to run the IAIPW estimator in its true form. The script 'sample_call.R' sources this folder and then runs the function. 

The asymptotic variance takes a long time to run. If it is not needed, you can remove the calls in the script 'IAIPW.R'. An option may be added at a later date. 
