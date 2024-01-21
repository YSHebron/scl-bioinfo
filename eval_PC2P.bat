set maindir=%CD%

echo.
echo ==========================================
echo ================ EVALUATING CLUSTERS =====
echo ==========================================
echo.
cd %maindir%\code\PC2P\EVALUATION
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
perl evaluate_clusters.pl -i "..\yeast_xval_output\clusters integrated"  -c "..\complexes_yeast.txt" -n 10 -l 0 -m 1 -x "..\xval_yeast.txt"

cd %maindir%