#  To run:
# .\pipeline.ps1 -i <inputPPIN> -r <goldStandard> -o <outputdirectory>
#  Example:
# .\pipeline.ps1 -i ".\data\Dummy_CYC_testonly.txt"  -r ".\data\Yeast\CYC_complexes.txt" -o ".\data\Results"
# .\pipeline.ps1 -i ".\data\Yeast\Collins\collins2007.txt"  -r ".\data\Yeast\CYC_complexes.txt" -o ".\data\Results\Dummy\Trial" -n ".\data\Negatome\negatome_2_mix_mapped.txt"
# --negfile data/Negatome/negatome_2_mix_mapped.txt

param (
    [Parameter(Mandatory=$true)] [string]$i,
    [Parameter(Mandatory=$true)] [string]$r,
    [Parameter(Mandatory=$true)] [string]$o,
    [string]$n
)

# Check if $i is an existing file
if (-not (Test-Path -Path $i -PathType Leaf)) {
    Write-Error "'$i' is not a valid file."
    exit 1
}

# Check if $r is an existing file
if (-not (Test-Path -Path $r -PathType Leaf)) {
    Write-Error "'$r' is not a valid file."
    exit 1
}

# Check if $o is an existing directory
if (-not (Test-Path -Path $o -PathType Container)) {
    Write-Output "The specified output directory '$o' does not exist. Creating it now."
    New-Item -Path $o -ItemType Directory
    if (-not (Test-Path -Path $o -PathType Container)) {
        Write-Error "Failed to create the output directory '$o'."
        exit 1
    }
}

Write-Output "Running P5COMP..."
Write-Output "PPIN: '$i'"
Write-Output "Ref:  '$r'"
Write-Output "Output:   '$o'"
Write-Output ""

# Developer parameters
$filteredfile="data/Interm/filtered_ppin.txt"
$decompfile="data/Interm/decomp_ppin.txt"
$hubfile="data/Interm/hub_proteins.txt"
$clusters_PC2P="data/Results/Dummy/PC2P_predicted.txt"
$iAdjustCD_outfile="data/Interm/ppin_adjusted.txt"

# Denoising -> data/Interm/filtered_ppin.txt
## Filtering: Negatome and PerProteinPair filtering.
## Note: This pipeline is packaged with Negatome 2.0 datasets.
python code/filtering.py $i $r $filteredfile --negfile $n --confidence 0.33

### DECOMP 1: Hub Removal
### -> data/Interm/ppin_adjusted.txt
### -> data/Interm/decomp_ppin.txt, data/Interm/hub_proteins.txt 
python code/iAdjustCD.py $filteredfile $iAdjustCD_outfile
python code/hub_remove.py $filteredfile $hubfile $decompfile

# Parallel Clustering
## 1. PC2P with hub return -> data/Results/Dummy/PC2P_predicted.txt, data/Results/Dummy/PC2P_postprocessed.txt
Write-Output "=========================================="
Write-Output "================= PC2P ==================="
Write-Output "=========================================="
Write-Output ""
$predictsfile_PC2P = Join-Path -Path $o -ChildPath "PC2P_predicted.txt"
$postprocessed_PC2P = Join-Path -Path $o -ChildPath "PC2P_postprocessed.txt"
python code/PC2P/PC2P.py $decompfile $predictsfile_PC2P -p mp
python code/hub_return.py $predictsfile_PC2P $iAdjustCD_outfile $hubfile $filteredfile $postprocessed_PC2P

Write-Output ===========================`n

## 2. CUBCO+
Write-Output "=========================================="
Write-Output "================ CUBCO+ =================="
Write-Output "=========================================="
Write-Output ""
$predictsfile_CUBCO = Join-Path -Path $o -ChildPath "CUBCO+_predicted.txt"
$postprocessed_CUBCO = Join-Path -Path $o -ChildPath "CUBCO+_postprocessed.txt"
python code/CUBCO+/CUBCO.py $decompfile $o $predictsfile_CUBCO
python code/hub_return.py $predictsfile_CUBCO $iAdjustCD_outfile $hubfile $filteredfile $postprocessed_CUBCO

Write-Output ===========================`n

## 3. ClusterOne
Write-Output "=========================================="
Write-Output "============== ClusterOne ================"
Write-Output "=========================================="
Write-Output ""
$jarPath = ".\code\ClusterOne\cluster_one-1.0.jar"
$outfile = "\ClusterOne_predicted.txt"
$predictsfile_ClusterOne = $o+$outfile
New-Item -Path $predictsfile_ClusterOne -ItemType File
$javaCommand = "java -jar ""$jarPath"" ""$filteredfile"" > ""$predictsfile_ClusterOne"""
Invoke-Expression $javaCommand
## Score clusters
$postprocessed_ClusterOne = Join-Path -Path $o -ChildPath "ClusterOne_postprocessed.txt"
python code/ClusterOne/cluster_one_scoring.py $filteredfile $predictsfile_ClusterOne $postprocessed_ClusterOne

Write-Output ===========================`n

# Ensemble Clustering
Write-Output "=========================================="
Write-Output "================ P5COMP =================="
Write-Output "=========================================="
Write-Output ""
$P5COMP_clusters = Join-Path -Path $o -ChildPath "P5COMP_clusters.txt"
python code/ensemble.py $postprocessed_ClusterOne $postprocessed_CUBCO $postprocessed_PC2P $P5COMP_clusters
Write-Output ===========================`n

# Evaluation
Write-Output "=========================================="
Write-Output "============== Evaluating ================"
Write-Output "=========================================="
Write-Output ""
python code/eval.py $postprocessed_ClusterOne $postprocessed_CUBCO $postprocessed_PC2P $P5COMP_clusters $r $o