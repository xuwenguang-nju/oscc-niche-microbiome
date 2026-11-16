# Convert a .docx to PDF with Word, in a background job so a modal Office dialog
# cannot block the caller forever.
#
# Why a job: an earlier attempt called Word synchronously and hung indefinitely.
# The cause is that Word, once killed with Stop-Process while it had an open document,
# shows a Document-Recovery / Safe-Mode dialog on the next start. That dialog is modal
# and never answers, so the synchronous COM call waits for ever. Running the COM work
# in a child job, with AutomationSecurity disabled and the document opened read-only,
# avoids the recovery prompt and lets us bound the wait.
param(
    [Parameter(Mandatory = $true)][string]$Docx,
    [Parameter(Mandatory = $true)][string]$Pdf,
    [int]$TimeoutSeconds = 90
)

$ErrorActionPreference = "Continue"

Get-Process WINWORD -ErrorAction SilentlyContinue | Stop-Process -Force
Start-Sleep -Seconds 2

$job = Start-Job -ArgumentList $Docx, $Pdf -ScriptBlock {
    param($d, $p)
    $w = New-Object -ComObject Word.Application
    $w.Visible = $false
    $w.DisplayAlerts = 0
    $w.AutomationSecurity = 3                      # msoAutomationSecurityForceDisable
    $doc = $w.Documents.Open($d, $false, $true)    # ConfirmConversions, ReadOnly
    $doc.SaveAs2($p, 17)                           # wdFormatPDF
    $doc.Close(0)
    $w.Quit()
    [System.Runtime.Interopservices.Marshal]::ReleaseComObject($w) | Out-Null
    "ok"
}

$done = Wait-Job $job -Timeout $TimeoutSeconds
if ($done) {
    Receive-Job $job | Write-Output
} else {
    "TIMEOUT after $TimeoutSeconds s"
}

Stop-Job $job -ErrorAction SilentlyContinue
Remove-Job $job -Force -ErrorAction SilentlyContinue
Get-Process WINWORD -ErrorAction SilentlyContinue | Stop-Process -Force

if (Test-Path $Pdf) {
    "pdf: {0:N0} KB" -f ((Get-Item $Pdf).Length / 1KB)
} else {
    "PDF NOT CREATED"
}
