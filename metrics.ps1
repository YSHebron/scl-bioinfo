# Set the output CSV file and write the header
$outputCsv = "results.csv"
"Method, GldStd, PPIN, Predicts, Refs, Precision, Recall, F1-score, F2-score, AUC-PR, MMR, Sensitivity, PPP, Accuracy, F-Match, Separation" | Out-File -FilePath $outputCsv -Encoding utf8

# Define the methods, gold standards, and PPINs
$methods = @("PC2P", "CUBCO+", "ClusterOne")

$gldstds = @("CYC", "SGD")

$ppins = @("Collins", "Gavin", "KroganCore", "KroganExt", "BIM")

# Initialize variables
$p = ""
$r = ""
$o = "data/Results/P5COMP"

foreach ($gldstd in $gldstds) {
    if ($gldstd -eq "CYC") {
        $r = "eval/CYC2008.txt"
    } elseif ($gldstd -eq "SGD") {
        $r = "eval/SGD.txt"
    }

    foreach ($ppin in $ppins) {
        switch ($ppin) {
            "Collins" { $p = "eval/Collins.txt" }
            "Gavin" { $p = "eval/Gavin.txt" }
            "KroganCore" { $p = "eval/KroganCore.txt" }
            "KroganExt" { $p = "eval/KroganExt.txt" }
            "BIM" { $p = "eval/BIM.txt" }
        }

        # PC2P, CUBCO+, and ClusterOne
        foreach ($method in $methods) {
            $attribs = "${method}-${gldstd}-${ppin}"
            & ./pipeline3.ps1 -ppinfile $p -reffile $r -outputdir $o -filter perpair -attribs $attribs

        }
    }
}
