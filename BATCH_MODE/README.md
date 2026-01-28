**This directory contains python and shell routines to execute multiple instances of CPMOD_mp in batch**. 

The main routine to be modified is run_batch.py which includes settings of the working directory and input settings relevant to contrails. The pbe solver executable is compiled for UNIX environment only, so run_batch.py cannot be executed from a Windows session, for example interactively using spyder. Please refer to the shell script to understand how to execute the python wrapper (run_batch.py) from a Linux terminal.

Currently, this version does not allow seamless integration with LES input data for the jet quantities. For this, the user has to modify the relevant subroutines in the src_mp/.F90 files and recompile cpmod_mp from compile_mp/Makefile
