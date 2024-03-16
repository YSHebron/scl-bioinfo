$code = ".\code\PC2P\PC2P.py"
$inputdir = ".\code\PC2P\Yeast\FilteredPPINs"
$outputdir = ".\code\PC2P\Results\FromFiltered"
$pool_thresh = 100

$ppins = (Get-ChildItem ".\code\PC2P\Yeast\FilteredPPINs" -Recurse).Name -replace '\..*$'

function inputfile($inputdir, $ppin, $gldstd) {
    ($inputdir + "\" + $ppin + "\" + $ppin + "_" + $gldstd + "_weighted.txt")
}

foreach ($ppin in $ppins) {
    python $code ($inputdir + "\" + $ppin + ".txt") $outputdir -p --pool_thresh $pool_thresh
    Write-Output ================
}

<#
$code = ".\code\PC2P\PC2P.py"
$inputdir = ".\code\PC2P\Yeast"
$outputdir = ".\code\PC2P\Results\"
$ppins = @('Collins','Gavin','KroganCore','KroganExt')
$gldstds = @('CYC', 'SGD')
$pool_thresh = 100

function inputfile($inputdir, $ppin, $gldstd) {
    ($inputdir + "\" + $ppin + "\" + $ppin + "_" + $gldstd + "_weighted.txt")
}

foreach ($ppin in $ppins) {
    foreach ($gldstd in $gldstds) {
        python $code (inputfile $inputdir $ppin $gldstd) $outputdir -p --pool_thresh $pool_thresh
        Write-Output ================
    }
}
#>