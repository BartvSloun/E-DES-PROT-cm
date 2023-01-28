S3 Appendix to the manuscript "Predicting the postprandial glucose response in tissue-specific insulin resistant phenotypes: a computational modeling approach".
E-DES-PROT model MATLAB implementation
Matlab version: R2018b

CONTENTS:
run_model.m 		- script to run model.

diffeq_diabetes.m 	- model equations
errorFunction.m 	- error function for optimization
integratorfunG.m 	- integrator function for integral part of PID insulin model
load_pHealthy_c.m 	- parameters and constants
load_x0_input.m		- initialize inputs
example_data.csv 	- example data set (artificial data) 

example_figure.png 	- running the run_model.m file should produce this plot. Note that slight differences due to default settings of the optimizer in different MATLAB versions might occur.

USE:
Open and execute the run_model.m file in Matlab.

For reproducibility, running run_model.m without changes should result in the figure example.png.