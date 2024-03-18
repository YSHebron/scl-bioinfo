param (
    [switch]$evalonly = $false
)

$code = ".\code\PC2P\PC2P.py"
$inputdir = ".\code\PC2P\Yeast\FilteredPPINs"
$outputdir = ".\code\PC2P\Results\FromFiltered"
$pool_thresh = 100

$ppins = (Get-ChildItem ".\code\PC2P\Yeast\FilteredPPINs" -Recurse).FullName

if (-Not $evalonly) {
    foreach ($ppin in $ppins) {
        python $code $ppin $outputdir -p --pool_thresh $pool_thresh
        Write-Output ===========================`n
    }
    
}

$predicteddirs = @("FilteredSequential", "FilteredMP", "FilteredRay")

foreach ($filtering in $predicteddirs) {
    $filestoeval = (Get-ChildItem ".\code\PC2P\Results\$filtering" -Recurse).FullName
    $filenamesonly = (Get-ChildItem ".\code\PC2P\Results\$filtering" -Recurse).Name
    Clear-Content ".\code\PC2P\Analysis\$filtering\auc_only.txt"
    for ($i=0; $i -lt $filestoeval.Length; $i++) {
        $gldstd = ($filenamesonly[$i] -split '_')[1]
        # python .\code\PC2P\PC2P_eval.py predictsfile complexfile outputdir
        python .\code\PC2P\PC2P_eval.py $filestoeval[$i] (".\code\PC2P\Yeast\$gldstd"+"_complexes.txt") ".\code\PC2P\Analysis\$filtering" | Out-File ".\code\PC2P\Analysis\$filtering\auc_only.txt" -Append
        Write-Output "" | Out-File ".\code\PC2P\Analysis\$filtering\auc_only.txt" -Append
    }
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