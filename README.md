# SBLF

1. C code to fit the spatial Bayesian latent factor (SBLF) model

2. RCode: R code for the analysis of model (SBLF) fitting results and model comparison with linear model and voxel-wise model.

The 'main2.c' files is for real data analysis, you can just run the code in vscode or on terminal.

# to compile the C code on terminal:
gcc -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas -lfftw3 -lm main2.c -o model2
./model2

# Note: 
1) please specify your own header and library search pathway if they are not found automatically. Here my header and library search pathway are '-I/usr/local/include' and '-L/usr/local/lib', respectively.
2) make sure that the GSL and FFTW libraries have been downloaded for use.

# Rcode:
'Functions.R': 
the functions are used for read data, fit linear and voxelwise regression models, comparing modeling results and draw figures.

'RealDataAnalysis.R': 
for real data analysis.

We use the a subset of randomly selected real data here. All of the imaging data is from ROI-4 (left amygdala region). The subset is split into a training (n=50) and a test (n=10) set.

Predictors: 
	'ROI_dat_mat.txt' (training) and 'ROI_dat_mat_test.txt' (test) 
	32 imaging predictors; 315 voxels per image.

Outcome:
	'ROI_task.txt' (training) and 'ROI_task_test.txt' (test)
	1 imaging outcome; 315 voxels per image.
	This outcome image is from the Emotion task domain. 
