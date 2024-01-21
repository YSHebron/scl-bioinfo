@echo off

set inputfile=%1%
shift
set outputdir=%1%
shift
set iters=%1%
shift
set mode=%1
shift
set pool_thresh=%1
shift
set num_procs=%1
shift

if "%mode%"=="" (
	echo Mode undefined
	exit /b
)

set maindir=%CD%

echo ==========================================
echo ================= PC2P ===================
echo ==========================================
echo. 

:: Example:
:: .\run_PC2P.ps1 .\code\PC2P\Human\PIPS\PIPS_Corum_Graph.txt .\code\PC2P\results\ 10 2
python code/PC2P/main.py %inputfile% %outputdir% %iters% %mode% %pool_thresh% %num_procs%

echo.
echo ==========================================
echo ================ EVALUATING CLUSTERS =====
echo ==========================================
echo.
cd %maindir%\EVALUATION
echo Evaluating PC2P...
:: Use -n 10 cluster iterations produced in "yeast_xval_output\clusters integrated"
:: inputdatafile = data_yeast.txt
:: complexfile = complexes_yeast.txt
:: goschemefile = go_scheme_all_mine.txt
:: goannotfile = my_go_associations.sgd.propagated.txt
:: xvalfile = xval_yeast.txt
:: outputdir = yeast_xval_output
:: paramset = yeast

:: -i input clusters file
:: -c complexes file
:: -n number of iterations
:: -x xval file, only when -n option is used
:: Other parameters: -m -l -k -s -a -u (see evaluate_clusters.pl)
perl evaluate_clusters.pl -i "..\%outputdir%\clusters integrated"  -c "..\%complexfile%" -n 10 -l 0 -m 1 -x "..\%xvalfile%" > "..\%outputdir%\results clusters integrated.txt"

cd %maindir%