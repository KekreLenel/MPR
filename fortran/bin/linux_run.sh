#!/bin/sh

debug_mode=0        # run fast (0) or debug (1)
make_prm_files=1    # create parameter files (1)
compile_exe=1       # compile executable (1)
solve_model=1       # solve model 
create_results=1    # run matlab results script (1) 

src_folder="../src/fortran/"

ulimit -s unlimited

# set variables for for OMP and NAG
export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=1
export NAG_KUSARI_FILE=lic.txt

# set environmental variables for intel fortran compiler
echo  "Set Intel Fortan Compiler vars."
source /opt/intel/oneapi/setvars.sh intel64 || exit 1
echo  "DONE"
echo  ""

# set environmental variables for nag library
nagfdir=/usr/local/NAG/nll6i271bl
echo  "Set NAG vars."
source ${nagfdir}/scripts/nagvars.sh -quiet int64 vendor dynamic || exit 1
echo  "DONE"
echo  ""

if [ $debug_mode -eq 0 ]
then
    # fast compilation 
    fcompile="ifort -O3 -xhost -qopenmp ${NAGLIB_FFLAGS} ${NAGLIB_FINCLUDE}"
else 
    # debug compilation 
    fcompile="ifort -O0 -init=snan -debug -traceback -check all,noarg_temp_created -nowarn -CB ${NAGLIB_FFLAGS} ${NAGLIB_FINCLUDE}"
fi 

# create output folder  
mkdir ../output
mkdir ../output/tmp

if [ $make_prm_files -eq 1 ]
then
    echo  "Create parameter files."
    rm ../src/params/*.csv
    matlab -nodisplay -nosplash -nodesktop <../src/params/create_param_files.m > ../output/tmp/param_output.txt 2> ../output/tmp/param_error.txt || exit 1
    echo  "DONE."
    echo  ""
fi 


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
    ${fcompile} ${src_folder}/base_lib.f90  ${src_folder}/mod_smolyak.f90 \
                ${src_folder}/mod_param.f90 ${src_folder}/mod_results.f90 \
                ${src_folder}/mod_calc.f90 ${src_folder}/mod_decomp.f90   \
                ${src_folder}/main.f90 -o main.exe ${NAGLIB_FLINK} || exit 1

    echo  "DONE."
    echo ""
fi 

# get number of calibrations to run
file="../output/tmp/n_comp.txt"
n_comp=$(cat "$file")
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

echo "Remove compilation files."
rm *.exe
rm *.o
rm *.mod
rm *_genmod.f90
echo  "DONE."
echo ""

if [ $create_results -eq 1 ]
then

    rm -r ../output/figures
    mkdir ../output/figures
    rm -r ../output/tables
    mkdir ../output/tables

    cp ../src/params/create_param_files.m ../output/tables 
     
    echo "Run MATLAB results script."
    matlab -nodisplay -nosplash -nodesktop <../src/matlab/main.m > ../output/tmp/results_output.txt 2> ../output/tmp/results_error.txt
    echo "DONE"

    pdflatex -shell-escape -output-directory=../output -halt-on-error ../output/results.tex | grep '^!.*' -A200 --color=always
    rm ../output/*.aux
    rm ../output/*.fdb_latexmk
    rm ../output/*.fls
    rm ../output/*.out
    rm ../output/*.log
fi


