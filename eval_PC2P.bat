@echo off

set maindir=%CD%

echo.
echo ==========================================
echo ================ EVALUATING CLUSTERS =====
echo ==========================================
echo.
echo Evaluating PC2P...

:: ./eval_PC2P.bat 
python code/PC2P/eval_PC2P.py code/PC2P/results/KroganExt/G_PredictedComplexes_iter0.txt

cd %maindir%