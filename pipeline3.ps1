# PowerShell script for running the pipeline

param (
    [string]$ppinfile,
    [string]$reffile,
    [string]$outputdir,
    [string]$negfile,
    [string]$filtering,
    [string]$attribs,
    [switch]$help
)

function Show-Help {
    Write-Output "usage: .\pipeline.ps1 [-p [ppinfile]] [-r [reffile]] [-o [outputdir]] [-n [negfile]] [-f [filter]] [-h]"
    Write-Output ""
    Write-Output "Runs P5COMP on the given PPIN file (ppinfile) and evaluates against the given gold standard (reffile)."
    Write-Output "Final predicted clusters will be written in outputdir."
    Write-Output "Important: Protein names (PID) should be in gene name (ordered locus) or KEGG format (ex. YLR075W) to match gold standards."
    Write-Output ""
    Write-Output "options:"
    Write-Output "    -p [ppinfile]       path to PPIN file (.txt) where each row is (u v w) (required)"
    Write-Output "    -r [reffile]        path to gold standard or reference complexes file (.txt) (required)"
    Write-Output "    -o [outputdir]      path to output directory (required)"
    Write-Output "    -n [negfile]        path to negatome (.txt) where each row is (u v) (optional)"
    Write-Output "    -f [filter]         filtering type (perpair or perprotein)"
    Write-Output "    -a [attribs]        attributes for evaluation file name of format 'algo-goldstd-ppin', ex: P5COMP-CYC-Collins"
    Write-Output "    -h                  show this help information"
}

function Log-Message {
    param (
        [string]$Message
    )
    Add-Content -Path "./debug.txt" -Value $Message
}

function Validate-File {
    param (
        [string]$FilePath
    )
    if (-Not (Test-Path -Path $FilePath -PathType Leaf)) {
        Write-Output "$FilePath is not a valid file."
        exit 1
    }
}

function Validate-Directory {
    param (
        [string]$DirPath
    )
    if (-Not (Test-Path -Path $DirPath -PathType Container)) {
        Write-Output "$DirPath is not a valid directory."
        exit 1
    }
}

if ($help) {
    Show-Help
    exit 0
}

if (-Not $ppinfile -or -Not $reffile -or -Not $outputdir) {
    Write-Output "Error: Missing -p, -r, and/or -o arguments. See help (-h)."
    exit 1
}

# Developer parameters
$filteredfile = "data/Interm/filtered_ppin.txt"

# Validate paths
Validate-File -FilePath $ppinfile
Validate-File -FilePath $reffile
Validate-Directory -DirPath $outputdir

# Create the output directory if it doesn't exist
if (-Not (Test-Path -Path $outputdir -PathType Container)) {
    New-Item -ItemType Directory -Path $outputdir | Out-Null
}

Write-Output "RUNNING INDEPENDENT CLUSTERING ALGOS:"
Write-Output "Denoising..."
python code/filtering.py $ppinfile $reffile $filteredfile --filtering $filtering

Write-Output "IND 1: Running PC2P..."
$predictsfile_PC2P = "$outputdir/PC2P_predicted.txt"
python code/PC2P/PC2P.py $filteredfile $predictsfile_PC2P -p mp

Write-Output "IND 2: Running CUBCO+..."
$predictsfile_CUBCO = "$outputdir/CUBCO+_predicted.txt"
python code/CUBCO+/CUBCO.py $filteredfile $outputdir $predictsfile_CUBCO

Write-Output "IND 3: Running ClusterOne..."
$predictsfile_ClusterOne = "$outputdir/ClusterOne_predicted.txt"
$postprocessed_ClusterOne = "$outputdir/ClusterOne_postprocessed.txt"
$jarPath = "code/ClusterOne/cluster_one-1.0.jar"
java -jar $jarPath $filteredfile > $predictsfile_ClusterOne

# Score ClusterOne clusters
python code/ClusterOne/cluster_one_scoring.py $ppinfile $predictsfile_ClusterOne $postprocessed_ClusterOne

# Evaluation (currently assumes running from root)
$PC2P_attribs = $x.Replace('-', '&')
python code/eval2.py $predictsfile_PC2P $reffile results.csv auc_pts.csv --attribs $attribs
python code/eval2.py $predictsfile_CUBCO $reffile results.csv auc_pts.csv --attribs $attribs
python code/eval2.py $postprocessed_ClusterOne $reffile results.csv auc_pts.csv --attribs $attribs
