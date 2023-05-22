#!/bin/bash
#SBATCH --job-name=KL
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=mlenel@princeton.edu
#SBATCH --output=KL_results.out
#SBATCH --error=KL_results.err
#SBATCH --mem=32000
module load matlab/R2019b

echo "Remove figures and tables"
rm -r ../output/figures
mkdir ../output/figures
rm -r ../output/tables
mkdir ../output/tables
echo "DONE"
cp ../src/params/create_param_files.m ../output/tables 
echo ""
echo "Create result files."
matlab -nodisplay -nosplash -nodesktop <../src/matlab/main.m > ../output/tmp/matlab_output_results.txt 2> ../output/tmp/matlab_error_results.txt
echo "DONE."


