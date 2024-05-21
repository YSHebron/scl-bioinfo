#  To run:
# .\pipeline.ps1 -i <inputPPIN> -r <goldStandard> -o <outputdirectory>
#  Example:
# .\pipeline.ps1 -i ".\data\Dummy_CYC_testonly.txt"  -r ".\data\Yeast\CYC_complexes.txt" -o ".\data\Results"
# .\pipeline.ps1 -i ".\data\Yeast\FilteredPPINs\Collins_CYC_perpair_weighted.txt"  -r ".\data\Yeast\CYC_complexes.txt" -o ".\data\Results"

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

# Denoising -> data/Interm/filtered_ppin.txt
## Filtering: Negatome and PerProteinPair filtering.
## Note: This pipeline is packaged with Negatome 2.0 datasets.
python code/filtering.py $i $r $filteredfile --negfile $n --confidence 0.33

## DECOMP 1: Hub Removal -> data/Interm/decomp_ppin.txt, data/Interm/hub_proteins.txt
python code/hub_remove.py $filteredfile $decompfile $hubfile

# Parallel Clustering



## 1. PC2P
Write-Output "=========================================="
Write-Output "================= PC2P ==================="
Write-Output "=========================================="
Write-Output ""
python code/PC2P/PC2P.py $i $o -p mp
#python code/PC2P/PC2P.py $decompfile $clusters_PC2P -p mp
Write-Output ===========================`n

## 2. CUBCO+
Write-Output "=========================================="
Write-Output "================ CUBCO+ =================="
Write-Output "=========================================="
Write-Output ""
python code/CUBCO+/CUBCO.py $i $o
Write-Output ===========================`n

## 3. ClusterOne
Write-Output "=========================================="
Write-Output "============== ClusterOne ================"
Write-Output "=========================================="
Write-Output ""
$jarPath = ".\code\ClusterOne\cluster_one-1.0.jar"
$outfile = "\clusters_clusterone.txt"
$outdir = $o+$outfile
New-Item -Path $outdir -ItemType File
$javaCommand = "java -jar ""$jarPath"" ""$i"" > ""$outdir"""
Invoke-Expression $javaCommand
## Score clusters
python ".\code\ClusterOne\cluster_one_scoring.py" $i $outdir
Write-Output ===========================`n

# Ensemble Clustering

# Evaluation
#python code/PC2P/PC2P_eval.py $clusters_PC2P $reffile $outputdir