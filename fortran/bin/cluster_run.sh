#!/bin/sh
debug_mode=0        # run fast (0) or debug (1)
make_prm_files=1    # create parameter files (1)
compile_exe=1       # compile executable (1)
solve_model=1       # solve model 
create_results=1    # run matlab results script (1) 
start_idx=1         # with which run to start
CPUS_PER_TASK=4

echo "Load Intel and Matlab"
module load intel openmpi
module load matlab
ulimit -s unlimited
echo "DONE"
echo ""

src_folder="../src/fortran/"
nagfdir=~/NAG/fll6i26dcl
nagldir="${nagfdir}/lib"
nagmkldir=~/NAG/fll6i26dcl/mkl_intel64_11.3.3/lib
nagmoddir="${nagfdir}/nag_interface_blocks"
flink="${nagldir}/libnag_nag.a \
-Wl,--start-group ${nagmkldir}/libmkl_intel_lp64.a \
${nagmkldir}/libmkl_intel_thread.a \
${nagmkldir}/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm"

if [ $debug_mode -eq 0 ]
then
    # fast compilation 
     fcompile="ifort -O3 -qopenmp -assume bscc -heap-arrays  -I${nagmoddir}"
else 
    # debug compilation 
    fcompile="ifort -O0 -init=snan -qopenmp -debug -traceback -check all,noarg_temp_created -warn all -CB -I${nagmoddir}"
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

	${fcompile}  ${src_folder}/base_lib.f90    ${src_folder}/mod_smolyak.f90 \
			     ${src_folder}/mod_param.f90   ${src_folder}/mod_calc.f90    \
			     ${src_folder}/mod_results.f90 ${src_folder}/mod_decomp.f90  \
                 ${src_folder}/main.f90 -o main.exe $flink || exit 1
 
    echo  "DONE."
	echo ""
fi 


# get number of calibrations to run
file="../output/tmp/n_comp.txt"
n_comp=$(cat "$file")
if [ $solve_model -eq 1 ]
then

echo "Prepare result folders."
for i in `seq $start_idx $n_comp`
do

    foo="../output/tmp/res_${i}"
    rm -r $foo
    mkdir $foo

echo "DONE"
echo ""
done
fi

echo "Submit job."
jid1=$(sbatch --array=1 -c $CPUS_PER_TASK --export=OMP_NUM_THREADS=$CPUS_PER_TASK --parsable ./job.slurm)
jid1b=$(sbatch --array=2-$n_comp --depend=afterany:$jid1 -c $CPUS_PER_TASK --export=OMP_NUM_THREADS=$CPUS_PER_TASK --parsable ./job.slurm)
echo "DONE."
echo ""


if [ $create_results -eq 1 ]
then
     
    echo "Run MATLAB results script."
    jid2=$(sbatch --time=1:30:00 --depend=afterany:$jid1b ./results.slurm)
    echo "DONE"

fi
