# Define arrays and values corresponding to metrics.sh
$methods = @("PC2P", "CUBCO+", "ClusterOne")
$gldstds = @("CYC", "SGD")
$ppins = @("Collins", "Gavin", "KroganCore", "KroganExt", "BIM")
$o = "data/Results/P5COMP"

# Initialize results.csv with headers
"Method, GldStd, PPIN, Predicts, Refs, Precision, Recall, F1-score, F2-score, AUC-PR, MMR, Sensitivity, PPP, Accuracy, F-Match, Separation" | Out-File -FilePath results.csv -Encoding utf8

foreach ($gldstd in $gldstds) {
    switch ($gldstd) {
        "CYC" { $r = "eval/CYC2008.txt" }
        "SGD" { $r = "eval/SGD.txt" }
    }
    foreach ($ppin in $ppins) {
        switch ($ppin) {
            "Collins" { $p = "eval/Collins.txt" }
            "Gavin" { $p = "eval/Gavin.txt" }
            "KroganCore" { $p = "eval/KroganCore.txt" }
            "KroganExt" { $p = "eval/KroganExt.txt" }
            "BIM" { $p = "eval/BIM.txt" }
        }
        # P5COMP
        .\pipeline2.ps1 -ppinfile $p -reffile $r -outputdir $o `
            -negfile data/Negatome/negatome_2_mix_mapped.txt -filtering perpair -attribs "P5COMP-${gldstd}-${ppin}"
        # # PC2P, CUBCO+, and ClusterOne
        # foreach ($method in $methods) {
        #     .\pipeline3.ps1 -ppinfile $p -reffile $r -outputdir $o -filtering perpair -attribs "${method}-${gldstd}-${ppin}"
        # }
    }
}
