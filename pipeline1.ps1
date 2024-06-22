# P5COMP: Parameter-free Pipeline for Predicting Problematic Protein Complexes
# Can only currently be run from repository root.
# Disclaimer: This still contain necessary usage parameters such as paths to files and directories
# Note that usage parameters are not clustering parameters.

# Sample: .\pipeline2.ps1 -p data\Yeast\Collins\collins2007.txt -r data\Yeast\CYC_complexes.txt -o data\Results\P5COMP -n data\Negatome\negatome_2_mix_mapped.txt -f perpair

param(
    [string]$ppinfile,
    [string]$reffile,
    [string]$outputdir,
    [string]$negfile,
    [string]$filtering,
    [string]$attribs
)

function Help {
@"
usage: .\pipeline2.ps1 [-p [ppinfile]] [-r [reffile]] [-o [outputdir]] [-n [negfile]] [-f [filter]] [-a [attribs]] [-h]

Runs P5COMP on the given PPIN file (ppinfile) and evaluates against the given gold standard (reffile).
Final predicted clusters will be written in outputdir.
Important: Protein names (PID) should be in gene name (ordered locus) or KEGG format (ex. YLR075W) to match gold standards.

options:
    -p [ppinfile]       path to PPIN file (.txt) where each row is (u v w) (required)
    -r [reffile]        path to gold standard or reference complexes file (.txt) (required)
    -o [outputdir]      path to output directory (required)
    -n [negfile]        path to negatome (.txt) where each row is (u v) (optional)
    -f [filter]         filtering type (perpair or perprotein)
    -a [attribs]        attributes for evaluation file name of format 'algo-goldstd-ppin', ex: P5COMP-CYC-Collins
    -h                  show this help information
"@
}

function Validate-File {
    param(
        [string]$file
    )
    if (-not (Test-Path $file -PathType Leaf)) {
        Write-Error "$file is not a valid file."
        exit 1
    }
}

function Validate-Directory {
    param(
        [string]$dir
    )
    if (-not (Test-Path $dir -PathType Container)) {
        Write-Error "$dir is not a valid directory."
        exit 1
    }
}

# if ($args.Count -eq 0) {
#     Write-Error "Error: No options supplied. See help below."
#     Help
#     exit 1
# }

if (-not $ppinfile -or -not $reffile -or -not $outputdir) {
    Write-Error "Error: Missing -p, -r, and/or -o arguments. See help (-h)."
    exit 1
}

# Create output directory if it doesn't exist
New-Item -ItemType Directory -Force -Path $outputdir | Out-Null

# Validate input files and directories
Validate-File $ppinfile
Validate-File $reffile
Validate-Directory $outputdir

Write-Output "Running P5COMP..."
Write-Output "PPIN: $(Resolve-Path $ppinfile)"
Write-Output "Ref: $(Resolve-Path $reffile)"
Write-Output "Output: $(Resolve-Path $outputdir)"

# Developer parameters
$filteredfile = "data/Interm/filtered_ppin.txt"
$decompfile = "data/Interm/decomp_ppin.txt"
$hubfile = "data/Interm/hub_proteins.txt"
$iAdjustCD_outfile = "data/Interm/ppin_adjusted.txt"

# Denoising -> data/Interm/filtered_ppin.txt
## Filtering: Negatome and either PerProteinPair / PerProtein filtering.
## Note: This pipeline is packaged with Negatome 2.0 datasets.
Write-Output "1. Denoising..."
python code/filtering.py $ppinfile $reffile $filteredfile --negfile $negfile --filtering $filtering

### DECOMP 1: Hub Removal
### -> data/Interm/ppin_adjusted.txt
### -> data/Interm/decomp_ppin.txt, data/Interm/hub_proteins.txt 
Write-Output "2. Hub Removal..."
python code/iAdjustCD.py $filteredfile $iAdjustCD_outfile
python code/hub_remove.py $filteredfile $hubfile $decompfile
# At this point, hubfile contains the hub proteins while decompfile contains the decomposed PPIN.
# The iAdjustCD_outfile contains the rescored PPIN, for use in hub_return.

Write-Output "Clustering Algorithms..."

# Parallel Clustering (Ensemble Clustering)
## 1. PC2P with hub return -> data/Results/Dummy/PC2P_predicted.txt, data/Results/Dummy/PC2P_postprocessed.txt
Write-Output "Running PC2P..."
$predictsfile_PC2P = Join-Path $outputdir "PC2P_predicted.txt"
$postprocessed_PC2P = Join-Path $outputdir "PC2P_postprocessed.txt"
python code/PC2P/PC2P.py $decompfile $predictsfile_PC2P -p mp
python code/hub_return.py $predictsfile_PC2P $iAdjustCD_outfile $hubfile $filteredfile $postprocessed_PC2P

## 2. CUBCO+ with hub return -> data/Results/Dummy/CUBCO+_predicted.txt, data/Results/Dummy/CUBCO+_postprocessed.txt
## Note: omit '+' character from varnames
Write-Output "Running CUBCO+..."
$predictsfile_CUBCO = Join-Path $outputdir "CUBCO+_predicted.txt"
$postprocessed_CUBCO = Join-Path $outputdir "CUBCO+_postprocessed.txt"
python code/CUBCO+/CUBCO.py $decompfile $outputdir $predictsfile_CUBCO
python code/hub_return.py $predictsfile_CUBCO $iAdjustCD_outfile $hubfile $filteredfile $postprocessed_CUBCO

## 3. ClusterOne with hub return -> data/Results/Dummy/ClusterOne_predicted.txt, data/Results/Dummy/ClusterOne_postprocessed.txt
Write-Output "Running ClusterOne..."
$predictsfile_ClusterOne = Join-Path $outputdir "ClusterOne_predicted.txt"
$postprocessed_ClusterOne = Join-Path $outputdir "ClusterOne_postprocessed.txt"
$jarPath = "code/ClusterOne/cluster_one-1.0.jar"
java -jar $jarPath $filteredfile > $predictsfile_ClusterOne

## Score clusters
python code/ClusterOne/cluster_one_scoring.py $ppinfile $predictsfile_ClusterOne $postprocessed_ClusterOne

# Ensemble Clustering
Write-Output "Finale: Running Ensemble Clustering..."
$final_clusters = Join-Path $outputdir "${attribs}_clusters.txt"
python code/ensemble.py $postprocessed_ClusterOne $postprocessed_CUBCO $postprocessed_PC2P $final_clusters

# Evaluation (currently assumes running from root)
python code/eval2.py $final_clusters $reffile results.csv auc_pts.csv --attribs $attribs
