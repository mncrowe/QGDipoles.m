This directory contains Matlab scripts for solving the N-layer quasi-geostrophic modon problem from `Modon solutions in an N-layer quasi-geostrophic model' by M. N. Crowe and E. R. Johnson.

Scripts for the three examples considered in this paper are provided, as is a general script (N_Layer_Modon.m) which can be modified to find any required solutions. These files have been tested on Matlab 2022b-2023b and are expected to work on some older versions. All scripts are written by M. N. Crowe based on theory devised by M. N. Crowe and E. R. Johnson. For comments, questions or suggestions, please email Matthew.Crowe2@newcastle.ac.uk.

Requirements:

 - Matlab (tested on 2022b-2023b but expected to work on most recent versions)
 - Matlab Optimization Toolbox

List of files:

functions/			- Directory with all functions required to determine modon solutions
	A_func.m		- Evaluates the function A_kj
	B_func.m		- Evaluates the function B_kj
	CalcPsi.m		- Calculates the streamfunctions and vorticities for a given solution
	create_domain.m		- Creates a 2D spatial domain with corresponding Fourier space domain
	DeltaN.m		- Defines operator used to find streamfunction from vorticitie, dependency of CalcPsi.m
	EVP.m			- Solves the one-parameter eigenvalue problem using matrix methods
	EVP_optim.m		- Solves the multi-parameter eigenvalue problem using root finding methods
	FT_2D.m			- Performs 2D Fourier transforms, dependency of CalcPsi.m
	jacobiP2.m		- Calculate Jacobi polynomials, faster that Matlab's `jacobiP', dependency of R_n.m
	JJ_int.m		- Evaluates a double Bessel function integral, dependency of A_func.m and B_func.m
	OrthogSpace.m		- Uses Gram-Schmidt to find basis of space orthogonal to inputs, dependency of EVP_optim.m
	R_n.m			- Evaluates the Zernike radial functions, dependency of CalcPsi.m
plotting/			- Directory with plotting scripts used in our examples
	cmap2.m			- Creates custom red/white/blue colorbar, dependency of Plot_2D.m
	dist_array.m		- Distributes input to multiple outputs, dependency of cmap2.m
	F_minmax.m		- Finds domain min and max of a given input, dependency of cmap2.m
	Plot_2D.m		- Plots a given 2D function with optional outline and circle centred on origin
N_Layer_Modon.m			- Example of an N-layer solution, currently set to N = 5
One_Layer_LRD.m			- One-layer LRD example in Section 4.1 of methods paper
Scaling_and_Accuracy_Tests.m	- Scaling and accuracy tests from supplementary material
Scaling_and_Accuracy_Tests_Par.m- Scaling and accuracy tests from supplementary material, parallelised
Three_Layer_Vortex.m		- Three-layer vortex example in Section 4.3 of methods paper
Two_Layer_Modon.m		- Two-layer modon example in Section 4.2 of methods paper
