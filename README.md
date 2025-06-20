CPMOD - Chemistry and Particulates MODelling
Version 1.0
11/05/2025

The code can be compiled with the makefile in directory compile. Options are offered for the gfortran compiler, either debug or release (optimised). Input files are in directories pbe and psr. Output files will be produced in directory pbe.
The code solves the population balance equation (PBE) in a spatially homogeneous environment, i.e. a perfectly stirred reactor (PFR).
To set up the problem, modify pbe.in for the PBE settings and psr.in for the PSR settings.
A number of test cases are provided, where the settings are built-in rather than read from input files. This ensures that they can always be used for reference if unaltered. Sample figures are provided in the directory test_cases_figures.

For more explanations on the methods in this code, you may refer to the following book and references therein:
Rigopoulos, S. ‘Population Balance of Particles in Flows: From Aerosols to Crystallisation’, Cambridge University Press, 2024.

This is the first version of the code.
The aggregation-fragmentation part can accommodate an arbitrary expanding grid, but a contracting grid is not yet implemented.
Various updates will be carried out in the future.

The following people contributed to the code:
Stelios Rigopoulos
Fabian Sewerin
Anxiong Liu
Binxuan Sun
Daniel O’Sullivan

Testing of the code:

Sample plots are provided for the test cases. To reproduce these plots, run the corresponding test case by selecting its option in the problem.in file. Output is in the following files in directory pbe:
psd.out - particle size distribution
psd_analytical.out - analytical solution for selected cases
psd_sp.out - self-preserving size distribution for Brownian aggregation

Instructions for plotting are shown below:

Options:
1. PSR with user inputs: currently set up to produce the same input as the test case of aggregation with the Brownian kernel (option 3).
2. Aggregation with the constant kernel: plot psd_analytical.out column 2 (y axis) vs column 1 (x axis) and psd.out column 3 (y axis) vs column 1 (x axis). This plots a comparison of the numerical and analytical solutions in terms of the number distribution, i.e. n(v) - v. See the sample plot agg_constant.pdf.
3. Aggregation with the Brownian kernel: plot psd_sp.out column 2 (y axis) vs column 1 (x axis). This plots the self-preserving distribution in terms of the similarity variable, i.e. psi - eta. See the sample plot agg_Brownian.pdf.
4. Aggregation-growth: plot psd_analytical.out column 3 (y axis) vs column 1 (x axis) and psd.out column 5 (y axis) vs column 1 (x axis). This plots a comparison of the numerical and analytical solutions in terms of the volume distribution, i.e. v*n(v) - v. See the sample plot agg-growth.pdf.
5. Nucleation-growth: plot psd.out column 3 (y axis) vs column 1 (x axis). See the sample plot nuc-growth.pdf.
6. Uniform binary fragmentation: plot psd_analytical.out column 2 (y axis) vs column 1 (x axis) and psd.out column 3 (y axis) vs column 1 (x axis). This plots a comparison of the numerical and analytical solutions in terms of the number distribution, i.e. n(v) - v. See the sample plot frac_un_bin.pdf.
7. Discrete PBE - aggregation with the Brownian kernel: plot psd_sp_discrete.out column 2 (y axis) vs column 1 (x axis). This plots the self-preserving distribution in terms of the similarity variable, i.e. psi - eta. See the sample plot agg_Brownian_discrete.pdf.



20/06/2025 - Luca Boscagli
This is a modified version of Stelios's book version that includes ice microphysics to model ice kinetics in a non-stationary, homogeneous environment
