@echo off
set inputdatafile=%1
shift
set complexfile=%1
shift
set goschemefile=%1
shift
set goannotfile=%1
shift
set xvalfile=%1
shift
set outputdir=%1
shift
set paramset=%1
shift
if "%paramset%"=="yeast" (
	set ngo=300
	set nhub=50
	set decompweight=.6
	set sssweight=1
) else if "%paramset%"=="human" (
	set ngo=1000
	set nhub=300
	set decompweight=.6
	set sssweight=.3
) else (
	set ngo=1000
	set nhub=300
	set decompweight=.6
	set sssweight=.3
)

set maindir=%CD%
echo input data = %inputdatafile%
echo complexes file = %complexfile%
echo GO scheme file = %goschemefile%
echo GO annotations file = %goannotfile%
echo X-val file = %xvalfile%
echo output dir = %outputdir%
echo Parameter set = %paramset%
echo ngo = %ngo%
echo nhub = %nhub%
echo DECOMP weight = %decompweight%
echo SSS weight = %sssweight%
echo.

if "%inputdatafile%"=="" (
	echo Input data file undefined
	exit /b
)
if "%complexfile%"=="" (
	echo Complex file undefined
	exit /b
)
if "%goschemefile%"=="" (
	echo GO scheme file undefined
	exit /b
)
if "%goannotfile%"=="" (
	echo GO annotations file undefined
	exit /b
)
if "%xvalfile%"=="" (
	echo X-val file undefined
	exit /b
)
if "%outputdir%"=="" (
	echo Output dir undefined
	exit /b
)
if "%ngo%"=="" (
	echo Ngo undefined
	exit /b
)
if "%nhub%"=="" (
	echo Nhub undefined
	exit /b
)
if "%decompweight%"=="" (
	echo DECOMP weight undefined
	exit /b
)
if "%sssweight%"=="" (
	echo SSS weight undefined
	exit /b
)

mkdir %outputdir%





echo ==========================================
echo ================== SWC ===================
echo ==========================================
echo. 
echo Scoring edges with SWC...
mkdir %outputdir%\SWC
cd %maindir%\SWC
perl score_edges.pl -i "..\%inputdatafile%" -c "..\%complexfile%" -m x -x "..\%xvalfile%" -a s -e 0 -l 2 -o "..\%outputdir%\SWC\swc" >> "%maindir%\%outputdir%\SWC\output swc.log"

echo.
echo Filtering SWC edges...
perl filter_scored_edges.pl -i "..\%outputdir%\SWC\swc scored_edges" -n 10 -s norm -k 20000 -o "..\%outputdir%\SWC\swc20k" >> "%maindir%\%outputdir%\SWC\output swc.log"



REM --------- Cluster SWC network --------
echo.
echo Extracting SWC clusters...

REM --- Cluster with CMC ---
echo.
echo Running CMC...
cd %maindir%\CLUSTERING\cmc
for /L %%i in (0,1,9) do cmc "..\..\%outputdir%\SWC\swc20k iter%%i.txt" 1 4 .5 .75 "..\..\%outputdir%\SWC\clusters_cmc swc20k iter%%i.txt" >> "%maindir%\%outputdir%\SWC\output swc.log"
for /L %%i in (0,1,9) do perl format_cmc_results.pl "..\..\%outputdir%\SWC\clusters_cmc swc20k iter%%i.txt" "..\..\%outputdir%\SWC\swc20k iter%%i.txt" 1 "..\..\%outputdir%\SWC\clusters_cmc_formatted swc20k iter%%i.txt" >> "%maindir%\%outputdir%\SWC\output swc.log"
for /L %%i in (0,1,9) do del "..\..\%outputdir%\SWC\clusters_cmc swc20k iter%%i.txt"

REM --- Cluster with ClusterOne ---
echo.
echo Running ClusterOne...
cd %maindir%\CLUSTERING\clusterone
for /L %%i in (0,1,9) do java -jar cluster_one-0.94.jar "..\..\%outputdir%\SWC\swc20k iter%%i.txt" -s 4 > "..\..\%outputdir%\SWC\clusters_clusterone swc20k iter%%i.txt"
for /L %%i in (0,1,9) do perl format_clusterone_results.pl "..\..\%outputdir%\SWC\clusters_clusterone swc20k iter%%i.txt" "..\..\%outputdir%\SWC\swc20k iter%%i.txt" 1 "..\..\%outputdir%\SWC\clusters_clusterone_formatted swc20k iter%%i.txt"  >> "%maindir%\%outputdir%\SWC\output swc.log"
for /L %%i in (0,1,9) do del "..\..\%outputdir%\SWC\clusters_clusterone swc20k iter%%i.txt"

REM --- Cluster with IPCA ---
echo.
echo Running IPCA...
cd %maindir%\CLUSTERING\ipca
for /L %%i in (0,1,9) do ConvertData.exe "..\..\%outputdir%\SWC\swc20k iter%%i.txt" PairScore Pair "..\..\%outputdir%\SWC\swc20k_ip%%i.txt"  >> "%maindir%\%outputdir%\SWC\output swc.log"
copy IPCA.exe "..\..\%outputdir%\SWC\"
cd %maindir%\%outputdir%\SWC\
for /L %%i in (0,1,9) do ipca -G"swc20k_ip%%i.txt" -S4 -P2 -T.6 -O"clus_i_swc20k_%%i.txt" >> "%maindir%\%outputdir%\SWC\output swc.log"
cd %maindir%\CLUSTERING\ipca
for /L %%i in (0,1,9) do perl format_ipca_results.pl "..\..\%outputdir%\SWC\clus_i_swc20k_%%i.txt" "..\..\%outputdir%\SWC\swc20k iter%%i.txt" 1 "..\..\%outputdir%\SWC\clusters_ipca_formatted swc20k iter%%i.txt" >> "%maindir%\%outputdir%\SWC\output swc.log"
del "..\..\%outputdir%\SWC\IPCA.exe"
for /L %%i in (0,1,9) do del "..\..\%outputdir%\SWC\swc20k_ip%%i.txt"
for /L %%i in (0,1,9) do del "..\..\%outputdir%\SWC\clus_i_swc20k_%%i.txt"

REM --- Cluster with Coach ---
echo.
echo Running Coach...
cd %maindir%\CLUSTERING\coach
for /L %%i in (0,1,9) do ConvertData.exe "..\..\%outputdir%\SWC\swc20k iter%%i.txt" PairScore Pair "..\..\%outputdir%\SWC\swc20k_coach iter%%i.txt" >> "%maindir%\%outputdir%\SWC\output swc.log"
cd %maindir%\%outputdir%\SWC
for /L %%i in (0,1,9) do "%maindir%\CLUSTERING\coach\coach" -in "swc20k_coach iter%%i.txt"  >> "%maindir%\%outputdir%\SWC\output swc.log"
cd %maindir%\CLUSTERING\coach
for /L %%i in (0,1,9) do perl format_coach_results.pl "..\..\%outputdir%\SWC\CoAch_complexes_swc20k_coach iter%%i.txt" "..\..\%outputdir%\SWC\swc20k iter%%i.txt" 1 "..\..\%outputdir%\SWC\clusters_coach_formatted swc20k iter%%i.txt"  >> "%maindir%\%outputdir%\SWC\output swc.log"
for /L %%i in (0,1,9) do del "..\..\%outputdir%\SWC\swc20k_coach iter%%i.txt"
for /L %%i in (0,1,9) do del "..\..\%outputdir%\SWC\CoAch_complexes_swc20k_coach iter%%i.txt"

REM --- Combine the clusters ---
echo.
echo Combining clusters...
cd %maindir%\CLUSTERING\combine
perl combine_cluster_results.pl -t .75 -n 10 -c "..\..\%outputdir%\SWC\clusters_cmc_formatted swc20k" -1 "..\..\%outputdir%\SWC\clusters_clusterone_formatted swc20k" -i "..\..\%outputdir%\SWC\clusters_ipca_formatted swc20k" -a "..\..\%outputdir%\SWC\clusters_coach_formatted swc20k" -o "..\..\%outputdir%\SWC\clusters_combined swc20k"  >> "%maindir%\%outputdir%\SWC\output swc.log"



echo.
echo ==========================================
echo ================== DECOMP ================
echo ==========================================
REM ------------- Extract top PPIs from data file ------
echo.
echo Extracting top PPIs...
cd %maindir%
mkdir %outputdir%\DECOMP
cd %maindir%\DECOMP
perl extract_data.pl -i "..\%inputdatafile%" -t "PPIREL" -n 20000 -o "..\%outputdir%\DECOMP\ppi20k.txt" >> "%maindir%\%outputdir%\DECOMP\output decomp.log"


REM ------------- Hub removal -------------
echo.
echo Performing hub removal...
cd %maindir%\DECOMP
perl remove_hubs.pl -i "..\%outputdir%\DECOMP\ppi20k.txt" -n %nhub% -h "..\%outputdir%\DECOMP\hubs.txt" -o "..\%outputdir%\DECOMP\PPIs_hub.txt" >> "%maindir%\%outputdir%\DECOMP\output decomp.log"

REM ------------- GO decomposition -------------
echo.
echo Performing GO decomposition...
cd %maindir%\DECOMP
perl get_decomp_goterms.pl -a "..\%goannotfile%" -s "..\%goschemefile%" -n %ngo% -o "..\%outputdir%\DECOMP\decompGOterms.txt" >> "%maindir%\%outputdir%\DECOMP\output decomp.log"
perl decomp_ppi_goterms.pl -i "..\%outputdir%\DECOMP\PPIs_hub.txt" -a "..\%goannotfile%" -d "..\%outputdir%\DECOMP\decompGOterms.txt" -o "..\%outputdir%\DECOMP\PPIs_hub_go" >> "%maindir%\%outputdir%\DECOMP\output decomp.log"

REM -- count number of GO decomposition terms --
findstr /R /N "^" "..\%outputdir%\DECOMP\decompGOterms.txt" | find /C ":" > "..\%outputdir%\DECOMP\numDecompTerms.tmp"
set /p numDecompTerms=<""..\%outputdir%\DECOMP\numDecompTerms.tmp"
del ""..\%outputdir%\DECOMP\numDecompTerms.tmp"
echo num GO Decomp terms = %numDecompTerms%
set /A highestDecompIndex=%numDecompTerms% - 1
echo highestDecompIndex = %highestDecompIndex%


REM ----- Cluster decomposed networks -----
echo.
echo Extracting DECOMP clusters...

REM --- Cluster with CMC ---
echo.
echo Running CMC...
cd %maindir%\CLUSTERING\cmc
for /L %%b in (0,1,%highestDecompIndex%) do cmc "..\..\%outputdir%\DECOMP\PPIs_hub_go_%%b.txt" 1 4 .5 .75 "..\..\%outputdir%\DECOMP\clusters_cmc PPIs_hub_go_%%b.txt" >> "%maindir%\%outputdir%\DECOMP\output decomp.log"
for /L %%b in (0,1,%highestDecompIndex%) do perl format_cmc_results.pl "..\..\%outputdir%\DECOMP\clusters_cmc PPIs_hub_go_%%b.txt" "..\..\%outputdir%\DECOMP\PPIs_hub_go_%%b.txt" 1 "..\..\%outputdir%\DECOMP\clusters_cmc_formatted PPIs_hub_go_%%b.txt" >> "%maindir%\%outputdir%\DECOMP\output decomp.log"
for /L %%b in (0,1,%highestDecompIndex%) do del "..\..\%outputdir%\DECOMP\clusters_cmc PPIs_hub_go_%%b.txt"
echo.
echo Recombining CMC's clusters from GO-decomposed subnetworks...
cd %maindir%\DECOMP
perl recombine_clusters.pl -i "..\%outputdir%\DECOMP\clusters_cmc_formatted PPIs_hub_go" -p "..\%outputdir%\DECOMP\ppi20k.txt" -n %numDecompTerms% -o "..\%outputdir%\DECOMP\clusters_cmc_recombined PPIs_hub_go.txt"
echo.
echo Re-adding hubs to CMC's clusters...
perl extract_data.pl -i "..\%inputdatafile%" -t "PPITOPO" -o "..\%outputdir%\DECOMP\ppitopo.txt" >> "%maindir%\%outputdir%\DECOMP\output decomp.log"
perl readd_hubs.pl -i "..\%outputdir%\DECOMP\clusters_cmc_recombined PPIs_hub_go.txt" -s "..\%outputdir%\DECOMP\ppitopo.txt" -h "..\%outputdir%\DECOMP\hubs.txt" -t 0.3 -o "..\%outputdir%\DECOMP\clusters_cmc_recombined_hubreadd PPIs_hub_go.txt" >> "%maindir%\%outputdir%\DECOMP\output decomp.log"




REM --- Cluster with ClusterOne ---
echo.
echo Running ClusterOne...
cd %maindir%\CLUSTERING\clusterone
for /L %%b in (0,1,%highestDecompIndex%) do java -jar cluster_one-0.94.jar "..\..\%outputdir%\DECOMP\PPIs_hub_go_%%b.txt" -s 4 > "..\..\%outputdir%\DECOMP\clusters_clusterone PPIs_hub_go_%%b.txt"
for /L %%b in (0,1,%highestDecompIndex%) do perl format_clusterone_results.pl "..\..\%outputdir%\DECOMP\clusters_clusterone PPIs_hub_go_%%b.txt" "..\..\%outputdir%\DECOMP\PPIs_hub_go_%%b.txt" 1 "..\..\%outputdir%\DECOMP\clusters_clusterone_formatted PPIs_hub_go_%%b.txt"  >> "%maindir%\%outputdir%\DECOMP\output decomp.log"
for /L %%b in (0,1,%highestDecompIndex%) do del "..\..\%outputdir%\DECOMP\clusters_clusterone PPIs_hub_go_%%b.txt"
echo.
echo Recombining ClusterOne's clusters from GO-decomposed subnetworks...
cd %maindir%\DECOMP
perl recombine_clusters.pl -i "..\%outputdir%\DECOMP\clusters_clusterone_formatted PPIs_hub_go" -p "..\%outputdir%\DECOMP\ppi20k.txt" -n %numDecompTerms% -o "..\%outputdir%\DECOMP\clusters_clusterone_recombined PPIs_hub_go.txt"
echo.
echo Re-adding hubs to ClusterOne's clusters...
perl readd_hubs.pl -i "..\%outputdir%\DECOMP\clusters_clusterone_recombined PPIs_hub_go.txt" -s "..\%outputdir%\DECOMP\ppitopo.txt" -h "..\%outputdir%\DECOMP\hubs.txt" -t 0.3 -o "..\%outputdir%\DECOMP\clusters_clusterone_recombined_hubreadd PPIs_hub_go.txt" >> "%maindir%\%outputdir%\DECOMP\output decomp.log"


REM --- Cluster with IPCA ---
echo.
echo Running IPCA...
cd %maindir%\CLUSTERING\ipca
for /L %%b in (0,1,%highestDecompIndex%) do ConvertData.exe "..\..\%outputdir%\DECOMP\PPIs_hub_go_%%b.txt" PairScore Pair "..\..\%outputdir%\DECOMP\ip_%%b.txt"  >> "%maindir%\%outputdir%\DECOMP\output decomp.log"
copy IPCA.exe "..\..\%outputdir%\DECOMP\"
cd "..\..\%outputdir%\DECOMP\"
for /L %%b in (0,1,%highestDecompIndex%) do ipca -G"ip_%%b.txt" -S4 -P2 -T.6 -O"clus_ipca_%%b.txt" >> "%maindir%\%outputdir%\DECOMP\output decomp.log"
cd %maindir%\CLUSTERING\ipca
for /L %%b in (0,1,%highestDecompIndex%) do perl format_ipca_results.pl "..\..\%outputdir%\DECOMP\clus_ipca_%%b.txt" "..\..\%outputdir%\DECOMP\PPIs_hub_go_%%b.txt" 1 "..\..\%outputdir%\DECOMP\clusters_ipca_formatted PPIs_hub_go_%%b.txt" >> "%maindir%\%outputdir%\DECOMP\output decomp.log"
del "..\..\%outputdir%\DECOMP\IPCA.exe"
for /L %%b in (0,1,%highestDecompIndex%) do del "..\..\%outputdir%\DECOMP\ip_%%b.txt"
for /L %%b in (0,1,%highestDecompIndex%) do del "..\..\%outputdir%\DECOMP\clus_ipca_%%b.txt"
echo.
echo Recombining IPCA's clusters from GO-decomposed subnetworks...
cd %maindir%\DECOMP
perl recombine_clusters.pl -i "..\%outputdir%\DECOMP\clusters_ipca_formatted PPIs_hub_go" -p "..\%outputdir%\DECOMP\ppi20k.txt" -n %numDecompTerms% -o "..\%outputdir%\DECOMP\clusters_ipca_recombined PPIs_hub_go.txt"
echo.
echo Re-adding hubs to IPCA's clusters...
perl readd_hubs.pl -i "..\%outputdir%\DECOMP\clusters_ipca_recombined PPIs_hub_go.txt" -s "..\%outputdir%\DECOMP\ppitopo.txt" -h "..\%outputdir%\DECOMP\hubs.txt" -t 0.3 -o "..\%outputdir%\DECOMP\clusters_ipca_recombined_hubreadd PPIs_hub_go.txt" >> "%maindir%\%outputdir%\DECOMP\output decomp.log"


REM --- Cluster with Coach ---
echo.
echo Running Coach...
cd %maindir%\CLUSTERING\coach
for /L %%b in (0,1,%highestDecompIndex%) do ConvertData.exe "..\..\%outputdir%\DECOMP\PPIs_hub_go_%%b.txt" PairScore Pair "..\..\%outputdir%\DECOMP\PPIs_hub_go_%%b_coach.txt" >> "%maindir%\%outputdir%\DECOMP\output decomp.log"
cd %maindir%\%outputdir%\DECOMP\
for /L %%b in (0,1,%highestDecompIndex%) do "%maindir%\CLUSTERING\coach\coach" -in "PPIs_hub_go_%%b_coach.txt"  >> "%maindir%\%outputdir%\DECOMP\output decomp.log"
cd %maindir%\CLUSTERING\coach
for /L %%b in (0,1,%highestDecompIndex%) do perl format_coach_results.pl "..\..\%outputdir%\DECOMP\CoAch_complexes_PPIs_hub_go_%%b_coach.txt" "..\..\%outputdir%\DECOMP\PPIs_hub_go_%%b.txt" 1 "..\..\%outputdir%\DECOMP\clusters_coach_formatted PPIs_hub_go_%%b.txt"  >> "%maindir%\%outputdir%\DECOMP\output decomp.log"
for /L %%b in (0,1,%highestDecompIndex%) do del "..\..\%outputdir%\DECOMP\PPIs_hub_go_%%b_coach.txt"
for /L %%b in (0,1,%highestDecompIndex%) do del "..\..\%outputdir%\DECOMP\CoAch_complexes_PPIs_hub_go_%%b_coach.txt"
echo.
echo Recombining Coach's clusters from GO-decomposed subnetworks...
cd %maindir%\DECOMP
perl recombine_clusters.pl -i "..\%outputdir%\DECOMP\clusters_coach_formatted PPIs_hub_go" -p "..\%outputdir%\DECOMP\ppi20k.txt" -n %numDecompTerms% -o "..\%outputdir%\DECOMP\clusters_coach_recombined PPIs_hub_go.txt"
echo.
echo Re-adding hubs to Coach's clusters...
perl readd_hubs.pl -i "..\%outputdir%\DECOMP\clusters_coach_recombined PPIs_hub_go.txt" -s "..\%outputdir%\DECOMP\ppitopo.txt" -h "..\%outputdir%\DECOMP\hubs.txt" -t 0.3 -o "..\%outputdir%\DECOMP\clusters_coach_recombined_hubreadd PPIs_hub_go.txt" >> "%maindir%\%outputdir%\DECOMP\output decomp.log"


REM --- Combine the clusters ---
echo.
echo Combining clusters...
cd %maindir%
cd CLUSTERING\combine
perl combine_cluster_results.pl -t .75 -n 1 -c "..\..\%outputdir%\DECOMP\clusters_cmc_recombined_hubreadd PPIs_hub_go.txt" -1 "..\..\%outputdir%\DECOMP\clusters_clusterone_recombined_hubreadd PPIs_hub_go.txt" -i "..\..\%outputdir%\DECOMP\clusters_ipca_recombined_hubreadd PPIs_hub_go.txt" -a "..\..\%outputdir%\DECOMP\clusters_coach_recombined_hubreadd PPIs_hub_go.txt" -o "..\..\%outputdir%\DECOMP\clusters_combined_recombined_hubreadd PPIs_hub_go.txt"  >> "%maindir%\%outputdir%\DECOMP\output decomp.log"


echo.
echo ==========================================
echo ================== SSS ===================
echo ==========================================
echo.
echo Scoring edges with SSS...
cd %maindir%
mkdir %outputdir%\SSS
cd %maindir%\SSS
perl score_edges_isoEM.pl -i "..\%inputdatafile%" -c "..\%complexfile%" -m x -x "..\%xvalfile%" -o "..\%outputdir%\SSS\sss" >> "..\%outputdir%\SSS\output sss.log"



echo Extracting small complexes...
perl extract_data.pl -i "..\%inputdatafile%" -t "PPIREL" -o "..\%outputdir%\SSS\ppi.txt" >> "%maindir%\%outputdir%\SSS\output sss.log"
for /L %%i in (0,1,9) do perl ExtractSmallComps.pl -i "..\%outputdir%\SSS\sss isoemiter2 scored_edges iter%%i.txt" -p "..\%outputdir%\SSS\ppi.txt" -r 10000 -o "..\%outputdir%\SSS\clusters_extract sss iter%%i.txt" >> "%maindir%\%outputdir%\SSS\output sss.log"



echo.
echo ==========================================
echo ================== INTEGRATE ALL ===========
echo ==========================================
echo.
echo Integrating generated clusters...
cd %maindir%
mkdir %outputdir%\INTEGRATE
cd %maindir%\INTEGRATE
for /L %%i in (0,1,9) do perl integrate_SWC_DECOMP_SSS.pl -1 "..\%outputdir%\SWC\clusters_combined swc20k iter%%i.txt" -2 "..\%outputdir%\DECOMP\clusters_combined_recombined_hubreadd PPIs_hub_go.txt" -3 "..\%outputdir%\SSS\clusters_extract sss iter%%i.txt" -y %decompweight% -z %sssweight% -t .75 -o "..\%outputdir%\clusters integrated iter%%i.txt" >> "%maindir%\%outputdir%\INTEGRATE\output integrate.log"




echo.
echo ==========================================
echo ================ EVALUATING CLUSTERS =====
echo ==========================================
echo.
cd %maindir%\EVALUATION
echo Evaluating combined approach (SWC + DECOMP + SSS)...
perl evaluate_clusters.pl -i "..\%outputdir%\clusters integrated"  -c "..\%complexfile%" -n 10 -l 0 -m 1 -x "..\%xvalfile%" > "..\%outputdir%\results clusters integrated.txt"

echo Evaluating SWC...
perl evaluate_clusters.pl -i "..\%outputdir%\SWC\clusters_combined swc20k"  -c "..\%complexfile%" -n 10 -l 0 -m 1 -x "..\%xvalfile%" > "..\%outputdir%\results clusters swc.txt"

echo Evaluating DECOMP...
for /L %%i in (0,1,9) do copy "..\%outputdir%\DECOMP\clusters_combined_recombined_hubreadd PPIs_hub_go.txt" "..\%outputdir%\DECOMP\clusters_combined_recombined_hubreadd PPIs_hub_go iter%%i.txt" 
perl evaluate_clusters.pl -i "..\%outputdir%\DECOMP\clusters_combined_recombined_hubreadd PPIs_hub_go"  -c "..\%complexfile%" -n 10 -l 0 -m 1 -x "..\%xvalfile%" > "..\%outputdir%\results clusters decomp.txt"
for /L %%i in (0,1,9) do del "..\%outputdir%\DECOMP\clusters_combined_recombined_hubreadd PPIs_hub_go iter%%i.txt" 

echo Evaluating SSS...
perl evaluate_clusters.pl -i "..\%outputdir%\SSS\clusters_extract sss"  -c "..\%complexfile%" -n 10 -l 0 -m 1 -x "..\%xvalfile%" > "..\%outputdir%\results clusters sss.txt"

cd %maindir%