#Load GFORTRAN for CPMOD
module load GCCcore

#Load python3
module load Python-bundle-PyPI

#Source python environment
source $HOME/python_cpmod/bin/activate

#Run CPMOD in batch mode
python3 run_batch.py
python3 run_postprocessing_batch.py
exit

