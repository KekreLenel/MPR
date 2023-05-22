@ECHO OFF 
set debug_mode=0
set make_prm_files=1    
set compile_exe=1       
set solve_model=1      
set create_results=1  
set src_folder=../src/fortran
set array_len = 10
@REM set OMP_NUM_THREADS=8
@REM set MKL_NUM_THREADS=1
@REM set NAG_KUSARI_FILE=lic.txt

echo  "Set Intel Fortan Compiler vars."
call "D:/Program Files (x86)/intel/setvars.bat" intel64 mod ilp64 || exit 1
echo  "DONE"
echo  ""
echo  "Set NAG vars."
call "D:\Program Files (x86)\NAG\NL27\nlw6i273el\batch\envvars.bat" || exit 1
echo  "DONE"
echo  ""

@REM if %debug_mode%==0 (set fcompile = "nagfor -abi=32 -ieee=full -compatible -I ") else (set fcompile="ifort -mkl -O0 -init=snan -debug -traceback -check all,noarg_temp_created -nowarn -CB %NAGLIB_FFLAGS% %NAGLIB_FINCLUDE%")

@REM if EXIST ..\output (rmdir /Q /S ..\output) else (echo "the file 'output' does not exist")
@REM echo  "removed"
@REM mkdir ..\output
@REM mkdir ..\output\tmp

@REM if %make_prm_files%==1 (
@REM     echo  "Create parameter files."
@REM     if EXIST ..\src\params\*.csv (del /Q /S ..\src\params\*.csv) else (echo "the param csv files does not exist")
@REM     matlab -nojvm -nosplash -nodesktop ../src/params/create_param_files || exit 1
@REM     echo  "DONE."
@REM     echo  ""
@REM ) 

if %compile_exe%==1 (
 	echo "Remove previous compilation files."
    @REM del /Q /S *.exe
    @REM del /Q /S .\*.o
    @REM del /Q /S .\*.mod
    @REM del /Q /S .\*_genmod.f90
    echo  "DONE."
 	echo ""

 	echo "Compile executable."

	ifort  /heap-arrays:50 %src_folder%/base_lib.f90 ^
			%src_folder%/mod_smolyak.f90 %src_folder%/mod_param.f90 ^
			%src_folder%/mod_calc.f90  %src_folder%/mod_results.f90 ^
			%src_folder%/mod_decomp.f90 %src_folder%/main.f90 ^
			nag_mkl_MT.lib mkl_rt.lib libiomp5md.lib user32.lib  -o main.exe   || exit 1
 

 	@REM ifort /iface:cvf %src_folder%/my.f90 -o my.exe || exit 1
 
    echo "DONE."
	echo ""
)

@REM # get number of calibrations to run
set file="../output/tmp/n_comp.txt"
set n_comp = 1
if %solve_model% == 1 ( 
    for %%i in (1, %n_comp%) do (
	rmdir /Q /S ..\output\tmp\res_%%i
	mkdir ..\output\tmp\res_%%i

 	echo "Run calibration %%i."
 	del /Q /S output.txt
    call main.exe %%i 
    echo  "DONE."
    )
)

echo "Remove compilation files."
@REM rm *.exe
@REM rm *.o
@REM rm *.mod
@REM rm *_genmod.f90

echo  "DONE."
echo ""

pause

