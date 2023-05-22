#!/bin/sh

debug_mode=0        # run fast (0) or debug (1)
make_prm_files=0    # create parameter files (1)
compile_exe=1       # compile executable (1)
solve_model=1       # solve model (1)
cleanup=0
create_results=0    # run matlab results script (1) 

src_folder="/Users/YourUsername/.../MPR/fortran/src/fortran"

ulimit -s 65532

# set variables for for OMP and NAG
export OMP_NUM_THREADS=3
export MKL_NUM_THREADS=1
export NAG_KUSARI_FILE=lic.txt

# set environmental variables for intel fortran compiler
echo  "Set Intel Fortan Compiler vars."
# source /opt/intel/bin/compilervars.sh intel64 || exit 1
source /opt/intel/oneapi/setvars.sh
echo  "DONE"
echo  ""

# set environmental variables for nag library
echo  "Set NAG vars."
nagfdir=/Users/YourUsername/NAG/nlmi627dbl
nagldir="${nagfdir}/ilp64/lib"
nagmkldir="${nagfdir}/mkl/lib"
nagmoddir="${nagfdir}lp64//nag_interface_blocks"
flink="${nagldir}/libnag_mkl.a \
       ${nagmkldir}/libmkl_intel_lp64.a \
       ${nagmkldir}/libmkl_intel_thread.a \
       ${nagmkldir}/libmkl_core.a \
       -framework IOKit -framework CoreFoundation \
       ${nagfdir}/rtl/lib/libiomp5.dylib -lpthread -lc++"
echo  "DONE"
echo  ""

# Compiling Command
if [ $debug_mode -eq 0 ]
then
    # fast compilation 
    fcompile="ifort -O3 -qopenmp -qmkl -I${nagmoddir} -init=snan"
else 
    # debug compilation 
    fcompile="ifort -O0  -qmkl -qopenmp -init=snan -debug -traceback -check all,noarg_temp_created -nowarn -CB -I${nagmoddir} -diag-disable 8889"
fi 

# create output folder  
mkdir ../output
mkdir ../output/tmp

# Create Parameter Files
if [ $make_prm_files -eq 1 ]
then
    echo  "Create parameter files."
    rm ../src/params/*.csv
    /Applications/MATLAB_R2021b.app/bin/matlab -nodisplay -nosplash -nodesktop <../src/params/create_param_files_New.m > ../output/tmp/param_output.txt 2> ../output/tmp/param_error.txt || exit 1
    echo  "DONE."
    echo  ""
fi 

# Compile
if [ $compile_exe -eq 1 ]
then
    echo "Remove previous compilation files."
    rm *.exe
    rm *.o
    rm *.mod
    rm *_genmod.f90
    echo  "DONE."
    echo ""

    echo "Compile executable."

    ${fcompile}  ${src_folder}/base_lib.f90 ${src_folder}/AuxCodes/mod_smolyak.f90  \
                 ${src_folder}/mod_param.f90 ${src_folder}/mod_calc.f90  \
                 ${src_folder}/mod_results.f90 ${src_folder}/mod_decomp.f90 \
                 ${src_folder}/main.f90 -o main.exe $flink || exit 1
 
    echo "Done"
 
    echo  "DONE."
    echo ""
fi 

# get number of calibrations to run
file="../output/tmp/n_comp.txt"
n_comp=$(cat "$file")

# Run Model
if [ $solve_model -eq 1 ]
then
# loop over calibrations
for i in `seq 1 $n_comp`
do
    # remove previous output files 
    foo="../output/tmp/res_${i}"
    rm -r $foo
    mkdir $foo

    echo "Run calibration ${i}."
    rm output.txt
    ./main.exe $i | tee ../output/tmp/output_${i}.txt
    echo  "DONE."
done
fi

if [ $cleanup -eq 1 ]
then

    echo "Remove compilation files."
    rm *.exe
    rm *.o
    rm *.mod
    rm *_genmod.f90
    rm -r *.exe.dSYM
    echo  "DONE."
    echo ""

fi

# Create Results
if [ $create_results -eq 1 ]
then

    rm -r ../output/figures
    mkdir ../output/figures
    rm -r ../output/tables
    mkdir ../output/tables

    cp ../src/params/create_param_files.m ../output/tables 
     
    echo "Run MATLAB results script."
    /Applications/MATLAB_R2021b.app/bin/matlab -nodisplay -nosplash -nodesktop <../src/matlab/result_scripts/create_results.m > ../output/tmp/results_output.txt 2> ../output/tmp/results_error.txt
    echo "DONE"


fi

