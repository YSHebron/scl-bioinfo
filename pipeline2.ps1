param(
    [string]$ppinfile,
    [string]$reffile,
    [string]$outputdir,
    [string]$filtering,
    [string]$attribs
)

function Help {
@"
usage: .\pipeline3.ps1 [-p [ppinfile]] [-r [reffile]] [-o [outputdir]] [-f [filter]] [-a [attribs]] [-h]

Runs P5COMP on the given PPIN file (ppinfile) and evaluates against the given gold standard (reffile).
Final predicted clusters will be written in outputdir.
Important: Protein names (PID) should be in gene name (ordered locus) or KEGG format (ex. YLR075W) to match gold standards.

options:
    -p [ppinfile]       path to PPIN file (.txt) where each row is (u v w) (required)
    -r [reffile]        path to gold standard or reference complexes file (.txt) (required)
    -o [outputdir]      path to output directory (required)
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

if ($args.Count -eq 0) {
    Write-Error "Error: No options supplied. See help below."
    Help
    exit 1
}

if (-not $ppinfile -or -not $reffile -or -not $outputdir) {
    Write-Error "Error: Missing -p, -r, and/or -o arguments. See help (-h)."
    exit 1
}

# Create output directory if it doesn't exist
New-Item -ItemType Directory -Force -Path $outputdir | Out-Null

Write-Output "RUNNING INDEPENDENT CLUSTERING ALGO:"
Write-Output "Denoising..."
python code/filtering.py $ppinfile $reffile $filteredfile --filtering $filtering

$method = $attribs.Split('-')[0]
switch ($method) {
    "PC2P" {
        Write-Output "IND 1: Running PC2P..."
        $predictsfile = Join-Path $outputdir "PC2P_predicted.txt"
        $postprocessed = Join-Path $outputdir "PC2P_postprocessed.txt"
        python code/PC2P/PC2P.py $filteredfile $predictsfile -p mp
        python code/PC2P/PC2P_scoring.py $filteredfile $predictsfile $postprocessed
        python code/eval2.py $postprocessed $reffile results.csv auc_pts.csv --attribs $attribs
    }
    "CUBCO+" {
        Write-Output "IND 2: Running CUBCO+..."
        $predictsfile = Join-Path $outputdir "CUBCO+_predicted.txt"
        $postprocessed = Join-Path $outputdir "CUBCO+_postprocessed.txt"
        python code/CUBCO+/CUBCO.py $filteredfile $outputdir $predictsfile
        python code/CUBCO+/cubco_scoring.py $filteredfile $predictsfile $postprocessed
        python code/eval2.py $postprocessed $reffile results.csv auc_pts.csv --attribs $attribs
    }
    "ClusterOne" {
        Write-Output "IND 3: Running ClusterOne..."
        $predictsfile = Join-Path $outputdir "ClusterOne_predicted.txt"
        $postprocessed = Join-Path $outputdir "ClusterOne_postprocessed.txt"
        $jarPath = "code/ClusterOne/cluster_one-1.0.jar"
        java -jar $jarPath $filteredfile > $predictsfile

        ## Score ClusterOne clusters
        python code/ClusterOne/cluster_one_scoring.py $ppinfile $predictsfile $postprocessed
        python code/eval2.py $postprocessed $reffile results.csv auc_pts.csv --attribs $attribs
    }
    default {
        Write-Error "Unknown method: $method"
        exit 1
    }
}
