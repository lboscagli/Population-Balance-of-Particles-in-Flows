## CPMOD - Chemistry and Particulates MODelling

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

## CPMOD FOR CONTRAILS 

20/06/2025 - Luca Boscagli

This is a modified version of Stelios's book version that includes ice microphysics to model ice kinetics in a non-stationary, homogeneous environment. Simplified depositional ice-crystal growth model based on [Karcher et al. (1996)](https://doi.org/10.1175/1520-0469(1996)053%3C3066:TICOJC%3E2.0.CO;2).


25/07/2025 - Luca Boscagli

Added functionality for condensational water droplet and depositional ice-crystal growth [Kärcher et al. (2015)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JD023491), and multi-particle (vPM and nvPM) activation [Ponsonby et al. 2025](https://doi.org/10.5194/egusphere-2025-1717). 

1. The new fortran routines I implemented are in the 'src_mp' directory and can be seamlessly compiled with the 'makefile' contained in the 'compile_mp' directory. 
2. The user can specify discretization of the PBE and characteristics of the nuclei of the particles via a user input file. This can be generated through the python script PBE_input_file.py contained in the 'pre_processing' directory. 
3. Results can be analysed via using the plotting functions contained in the python post-processing routines provided in the 'post_processing' directory.

### Repository structure
```text
│   cpmod_mp
│   LICENSE.txt
│   problem.in
│   README.md
│   run_cpmod.sh
│
├───BATCH_MODE
│       INPS_WRITER.py
│       PLOT_STATISTICS.py
│       README.md
│       run_batch.py
│       run_cpmod_batch.sh
│       run_postprocessing_batch.py
│
├───compile
│       agg_cfv_mod.mod
│       constants.mod
│       CPMOD_pbe.o
│       frag_cfv_mod.mod
│       gauss_mod.mod
│       ice_microphys_mod.mod
│       makefile
│       module_ice_microphysics.o
│       module_thermo.o
│       optimizer.mod
│       PBE_agg.o
│       PBE_agg_kernels.o
│       PBE_discrete.o
│       pbe_discrete_mod.mod
│       PBE_frag.o
│       PBE_general.o
│       PBE_growth.o
│       pbe_mod.mod
│       PBE_solver.o
│       PBE_test.o
│       PSR_pbe.o
│       thermo.mod
│
├───compile_mp
│       makefile
│
├───Info_material
│       combustion_water_vapor_formulas.docx
│       CPMOD_info.pdf
│       FUEL_water_vapor_concentration.py
│
├───pbe
│       pbe.in
│       pbe.in.default
│       pbe.in.nasa_buthanol_exp
│       psd.out
│       psd_analytical.out
│       psd_sp.out
│
├───post_processing
│       READ_MAT_psd_statistics.py
│       READ_psd.py
│       READ_psd_time_history.py
│
├───pre_processing
│       PBE_input_file_lognormal.py
│       PBE_input_file_monodisperse.py
│
├───psr
│       ice.in
│       ice.in.nasa_exp_buthanol
│       ice_nucleating_particles.in
│       psr.in
│       psr.in.nasa_exp_buthanol
│
├───src
│       CPMOD_pbe.f90
│       module_ice_microphysics.f90
│       module_thermo.f90
│       PBE_agg.f90
│       PBE_agg_kernels.f90
│       PBE_discrete.f90
│       PBE_frag.f90
│       PBE_general.f90
│       PBE_growth.f90
│       PBE_solver.f90
│       PBE_test.f90
│       PSR_pbe.f90
│
├───src_mp
│       CPMOD_pbe.f90
│       module_ice_microphysics.f90
│       module_thermo.f90
│       PBE_agg.f90
│       PBE_agg_kernels.f90
│       PBE_discrete.f90
│       PBE_frag.f90
│       PBE_general.f90
│       PBE_growth.f90
│       PBE_solver.f90
│       PBE_solver.modisperse.f90
│       PBE_test.f90
│       PSR_pbe.f90
│
└───test_cases_figures
        agg-growth.pdf
        agg_Brownian.pdf
        agg_Brownian_discrete.pdf
        agg_constant.pdf
        frag_un_bin.pdf
        nuc-growth.pdf
```text



