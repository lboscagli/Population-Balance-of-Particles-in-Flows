This directory contains python and shell routines to execute multiple instances of CPMOD_mp in batch. 

The main routine to be modified is run_batch.py which includes settings of the working directory and input settings relevant to contrails.

Currently, this version does not allow seamless integration with LES input data for the jet quantities. For this, the user has to modify the relevant subroutines in the src_mp/.F90 files and recompile cpmod_mp from compile_mp/Makefile
