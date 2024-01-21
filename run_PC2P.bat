@echo off

set inputfile=%1%
shift
set outputdir=%1%
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
:: ./run_PC2P.bat .\code\PC2P\Human\PIPS\PIPS_Corum_Graph.txt .\code\PC2P\results\ 2
python code/PC2P/main.py %inputfile% %outputdir% %mode% %pool_thresh% %num_procs%
