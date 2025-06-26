


First load the library "devtools" and hit load_all(). All functions and inputs will be loaded. Plots of the figures in the paper can be reproduced by going into the folder "Manuscript_Figrues".

Note that we have
not included the output of the simulation for space reasons, but have included RDS files of the kin-structures
produced from the simulation.

All code used is as follows:

The folder data/ contains files which contain UK age-specific fertility and mortality, downloaded from the HFD,HMD.


The folder R/ stores all of the theoretical model code:

1) age_mother_func.R infers the probable ages of ancestors in the population (see Section 2.3 in text)
2) fert_pdfs_Func.R transforms the fertility rates in a Leslie matrix into age-specific probability mass functions for offspring number
3) nth_convolution_Function.R and nth_fast_fourier_transform.R calculate the convolution of a sequence of distributions (as in Eq 3 in text)
4) Q_matrix_Func.R constructs the matrix F defined in Eq 12 in text (Section 2.5)
5) U_matrix_Func.R constructs the matrix U defined in Section 2.4
6) utility_funcs.R are some useful operations frequently used

The remainder of the folder -- ending in _PMF.R are functions which project age-specific probability mass functions of kin


The folder Sim/ contains code to run the simulation using Python's class attribute

1) BoidClass_PDF.py defines the Class from which we extract kinship in a natural population process
2) run.py runs the simulation

The folder Manuscript_Figures/ stores code to reproduce all figures in the manuscript: each file states which grahpic.

