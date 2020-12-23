#!/opt/perl5/bin/perl -w

$BIN     = '/home/wen/CorePhase/bin/';
$InstrFile = "../src.sac";
#$InstrFile = "/home/wen/Instr/LPWWSS.T\*1";


$ConvlvSource = "no";
$ConvlvInstr  = "yes";

@Recivers = (120, 121, 122, 123, 124, 125, 126, 127, 128, 129);

$RAD   = 111.194924748;

for($idist = 0; $idist < @Recivers; $idist++){
    #Reading out the Green's functions from library
    $dist  = $Recivers[$idist] * $RAD;
	#doing Cut and Parse 
	for($comp = 1; $comp < 2; $comp++){
	    if($comp == 0){
	        $sacfile   = "hy$Recivers[$idist].u";
	    } else{
	        $sacfile   = "hy$Recivers[$idist].w";
            }

	    $newfile   = "../$sacfile.prem";

            $sacfile = "/export/home/lehmann/wen/src/Dbounce/_p.sac$Recivers[$idist]";

            if($ConvlvSource eq "yes") {
		`$BIN/earqt if=$sacfile of=$newfile -S$DT1,$DT2,$DT3,1v`;
            }
            if($ConvlvInstr eq "yes") {
		system("$BIN/earqt if=$sacfile of=$newfile -I$InstrFile");
            }
        }
 }
    


