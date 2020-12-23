
#!/opt/perl5/bin/perl -w
###############################################################
#                                                             #
#   perl shell for running the hybrid method and Kirchhoff    #
#                                                             #
###############################################################

$BIN     = '/work/common/home/chl/wen/psvhy/my0/src/';
C$BIN     = '/work/common/home/chl/wen/psvhy/my0/src/';
$SACADD  = 'sacadd';
$PARFILE = "testiasp91.par";

$GenFDmodel           = "no";
$GRTrun               = "yes";
$FDrun                = "yes";
$KIRrun               = "no";
$JointGRTFDKIRrun     = "no";
$KIRGRTrun            = "no";
$WKMKIRrun            = "no";

# epicentral distance in degree
@Recivers = (130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145);

$WKMrun       = "no";
$WKMsolution  = "no";

$domewidth = 1500;
$low       = 10;
$center    = 500;
$heigh     = 250;
$xstart    = 10;
$xend      = 2440;
$zbase     = 258.2228; 
$vp        = 24.88374;
$vs        = 13.1866;
$d         = 10.07638;
$pfac      = 0.9;
$sfac      = 1.0;
$dfac      = 1.0;
$domevp    = $vp * $pfac;
$domevs    = $vs * $sfac;
$domed     = $d * $dfac;

$undulate1  = 0;
$undulate   = 0;
$undulate2  = 0;
$dn         = 2;
$np2        = 9;
$np         = 40;
$np1        = 12;
$fdmodel = &findpar($PARFILE, "fdmodel");


if($GenFDmodel eq "yes"){
    &GenDomeModel($fdmodel);
}

if($GRTrun eq "yes"){
    &GRTrun($PARFILE);
}

if($FDrun eq "yes"){
    &FDrun($PARFILE);
}

$RAD   = 111.194924748;
$break  = 0;
$npoint = 100;
$bstart = 0;
if($KIRrun eq "yes"){
     &KIRrun();
}

if($KIRGRTrun eq "yes"){
    &KIRGRTrun();
}


if($JointGRTFDKIRrun eq "yes"){
    &JointGRTFDKIRrun();
}
if($WKMKIRrun eq "yes"){
    &WKMKIRrun();
}


sub GenDomeModel{
    my($fdmodel) = @_;
    open(FDMODEL, ">$fdmodel") || die "Open FD model $fdmodel Failed\n";
    @l = `/home/wen/CorePhase/GenRuffGauss width=$domewidth heigh=$heigh xstart=$xstart zbase=$zbase xend=$xend vp=$vp vs=$vs d=$d domevp=$domevp domevs=$domevs domed=$domed undulate=$undulate undulate1=$undulate1 undulate2=$undulate2 dn=$dn np=$np np1=$np1 np2=$np2 low=$low center=$center`;
    print FDMODEL @l;
    close(FDMODEL);
}

sub GRTrun{
    my($par) = @_;
    #running the GRT calculations for input wavefields
    @l = `$BIN/aserpsvfd_my par=$par `;
    print @l;

    #demult the wavefield from GRT to input of FD
    system ("$BIN/demult par=$par  ");
}

sub FDrun{
    my($par) = @_;
    #running FD
    system ("$BIN/psvfd_my par=$par");
}

sub KIRrun{

    #demult the wavefield from FD to input of Kirchhoff 
    system ("$BIN/demult2kir par=$PARFILE");

    $nt_kir = &findpar($PARFILE, "nt_kir");
    $dt_kir = &findpar($PARFILE, "dt_WKM");
    $nx     = &findpar($PARFILE, "nx");
    $nh     = &findpar($PARFILE, "nh");
    $kdx    = &findpar($PARFILE, "kdx");
    $Nrecv  = int (($nx - $nh)/ $kdx / $npoint);
    $nt_kir *= 3;
    $greenfile_kir = &findpar($PARFILE, "greenfile_kir");

    for($idist = 0; $idist < @Recivers; $idist++){
        #Reading out the Green's functions from library
	$dist  = $Recivers[$idist] * $RAD;
	`echo dist=$dist`;
	system("$CBIN/Read_kirgreen par=$PARFILE dist=$dist");
	#doing Kirchhoff
	for($comp = 1; $comp < 2; $comp++){
	    if($comp == 0){
	        $fdkirfile = &findpar($PARFILE, "kirfile_x");
	        $sacfile   = "$Recivers[$idist].accx";
	    } else{
	        $fdkirfile = &findpar($PARFILE, "kirfile_z");
	        $sacfile   = "$Recivers[$idist].accz";
            }
	    $tstart = `$CBIN/kirch par=$PARFILE accfile=acc_junk dist=$dist fdkirfile=$fdkirfile breakup=$break nfre=$npoint breakupstart=$bstart`;
	    chomp($tstart);

	    system("$CBIN/line2point par=$PARFILE accfile=acc_junk nt_kir=$nt_kir");
	    `mv syn.SAC.1 $sacfile`;
	    &ChgSacHead($sacfile,"b",$tstart);

            if($break != 0){
	        system("$CBIN/line2point par=$PARFILE accfile=junkpars nrecv=$Nrecv nt_kir=$nt_kir");
	        for($k = 1; $k <= $Nrecv; $k++){
		     $file = "$sacfile.$k";
	            `mv syn.SAC.$k $file`;
	             &ChgSacHead($file,"b",$tstart);
	        }
	    }

            unlink "acc_junk";

        }
     }
     unlink  $greenfile_kir;
     unlink  $fdkirfile;
}

sub KIRGRTrun{
    $GREENSONLY = 0;
    $KIRONLY    = 1;
    $nt_kir    = &findpar($PARFILE, "nt_kir");
    $nt_kir *= 3;
    $dt_kir    = &findpar($PARFILE, "dt_WKM");
    $h         = &findpar($PARFILE, "h");
    $distmin   = &findpar($PARFILE, "distmin");
    $distmax   = &findpar($PARFILE, "distmax");
    $GRTtstart = &findpar($PARFILE, "GRTtstart");
    @xmin      = split(',',$distmin);
    @xmax      = split(',',$distmax);

    $greenfile_kir = &findpar($PARFILE, "greenfile_kir");

    for($i = 0; $i < @xmin; $i++){ 
	$xmin = $xmin[$i];
	$xmax = $xmax[$i];
	$nx   = int (($xmax - $xmin)/$h) + 3;
	$nh   = 0;
	$kdx  = 1;
        $Nrecv  = int (($nx - $nh)/ $kdx / $npoint) + 1;

	print "xmin=$xmin $xmax $nx $GRTtstart\n";

	#doing GRT
	for($comp = 0; $comp < 1; $comp++){
	    if($comp == 1){
	        $fdkirfile = &findpar($PARFILE, "GRTkirfile_x");
	        $sacfile   = "grt$Recivers[$idist].accx";
	    } else{
	        $fdkirfile = &findpar($PARFILE, "GRTkirfile_z");
	        $sacfile   = "grt$Recivers[$idist].accz";
            }
            if($KIRONLY != 1){
                system("$CBIN/aseries par=$PARFILE GRTkirfile=$fdkirfile comp=$comp");
            }
        }

        if($GREENSONLY != 1){
        for($idist = 0; $idist < @Recivers; $idist++){
            #Reading out the Green's functions from library
	    $dist  = $Recivers[$idist] * $RAD;
	    print "dist=$dist\n";
	    system("$CBIN/Read_kirgreen par=$PARFILE xmin=$xmin nx=$nx nl=0 dist=$dist nh=0");
	    #doing GRT and Kirchhoff 
	    for($comp = 0; $comp < 1; $comp++){
	        if($comp == 1){
	            $fdkirfile = &findpar($PARFILE, "GRTkirfile_x");
	            $sacfile   = "grt$Recivers[$idist].accx";
	        } else{
	            $fdkirfile = &findpar($PARFILE, "GRTkirfile_z");
	            $sacfile   = "grt$Recivers[$idist].accz";
                }

		#doing Kirchhoff 
	        $tstart = `$CBIN/kirch par=$PARFILE accfile=acc_junk dist=$dist xmin=$xmin nx=$nx nl=0 nh=$nh tstart=$GRTtstart fdkirfile=$fdkirfile breakup=$break nfre=$npoint breakupstart=$bstart`;
	        chomp($tstart);

	        system("$CBIN/line2point par=$PARFILE accfile=acc_junk nt_kir=$nt_kir");
	        `mv syn.SAC.1 $sacfile`;
	        &ChgSacHead($sacfile,"b",$tstart);

                if($break != 0){
	            system("$CBIN/line2point par=$PARFILE accfile=junkpars nrecv=$Nrecv nt_kir=$nt_kir");
	            for($k = 1; $k <= $Nrecv; $k++){
		         $file = "$sacfile.$k";
	                 `mv syn.SAC.$k $file`;
	                  &ChgSacHead($file,"b",$tstart);
	            }
	        }
            }
        }
	}
     }
     unlink $greenfile_kir;
}

sub ChgSacHead{
    my($file,$head,$val) = @_;
    open(SAC, "|sac2000");
    print SAC "r $file\n"; 
    print SAC "ch $head $val\n"; 
    print SAC "w over\n"; 
    print SAC "quit\n";
    close(SAC);
}


sub JointGRTFDKIRrun{

    $GRTOUT = '/home/wen/PKP/HY/PREM/';
    $DIR  = './';
    
    #segments are going to be sumed in the GRT calculation
    @segs = (1,2,3,4);

    for($idist = 0; $idist < @Recivers; $idist++){
        for($comp = 1; $comp < 2; $comp++){
	    if($comp == 0){
		$sac1  = "${GRTOUT}grt$Recivers[$idist].accx";
		$sac2  = "${DIR}$Recivers[$idist].accx";
		$sac3  = "${DIR}hy$Recivers[$idist].u";
            } else {
		$sac1  = "${GRTOUT}grt$Recivers[$idist].accz";
		$sac2  = "${DIR}$Recivers[$idist].accz";
		$sac3  = "${DIR}hy$Recivers[$idist].w";
            }
        }
	$in1 = $sac2;
	$d1  = 0.45; #this shift is due to the Gauss source	
	for($iseg = 0; $iseg < @segs; $iseg++){
	    $in2 = "$sac1.$segs[$iseg]";
	    `$SACADD if1=$in1 if2=$in2 of=$sac3 d1=$d1`;
	    $d1 = 0.0;
	    $in1 = $sac3;
        }
    }
}
sub WKMKIRrun{
    $nt_kir = &findpar($PARFILE, "nt_kir");
    $nt_kir *= 3;
    $dt_kir = &findpar($PARFILE, "dt_WKM");
    $nx     = &findpar($PARFILE, "nx");
    $nh     = &findpar($PARFILE, "nh");
    $kdx    = &findpar($PARFILE, "kdx");
    $Nrecv  = int (($nx - $nh)/ $kdx / $npoint);

    for($idist = 0; $idist < @Recivers; $idist++){
        #Reading out the Green's functions from library
	$dist  = $Recivers[$idist] * $RAD;
	`echo dist=$dist`;
	system("$CBIN/Read_WKMkirgreen par=$PARFILE greenfile_WKM=/home/u1/wen/Kirchhoff/test/GreenWKM RayParFile=/home/u1/wen/Kirchhoff/test/raypar dx_WKM=0.5 dist=$dist");
	#doing Kirchhoff
	for($comp = 0; $comp < 1; $comp++){
	    if($comp == 0){
	        $fdkirfile = &findpar($PARFILE, "kirfile_x");
	        $sacfile   = "WKM/$Recivers[$idist].accx";
	    } else{
	        $fdkirfile = &findpar($PARFILE, "kirfile_z");
	        $sacfile   = "WKM/$Recivers[$idist].accz";
            }
	    $tstart = `$CBIN/kirch par=$PARFILE accfile=acc_junk dist=$dist fdkirfile=$fdkirfile breakup=$break nfre=$npoint breakupstart=$bstart`;
	    chomp($tstart);

	    system("$CBIN/line2point par=$PARFILE accfile=acc_junk nt_kir=$nt_kir");
	    `mv syn.SAC.1 $sacfile`;

	    open(SAC, "|sac2000");
	    print SAC "r $sacfile\n"; 
	    print SAC "ch b $tstart\n"; 
	    print SAC "w over\n"; 
	    print SAC "quit\n";
	    close(SAC);

	    system("$CBIN/line2point par=$PARFILE accfile=junkpars nrecv=$Nrecv nt_kir=$nt_kir");
	    for($k = 1; $k <= $Nrecv; $k++){
	        `mv syn.SAC.$k $sacfile.$k`;

	        open(SAC, "|sac2000");
	        print SAC "r $sacfile.$k\n"; 
	        print SAC "ch b $tstart\n"; 
	        print SAC "w over\n"; 
	        print SAC "quit\n";
	        close(SAC);
	    }

        }
    }
}


###################################################################
#subroutine for finding the value of specific par parameter from  #
#par source file USAGE:                                           #
#$val = findpar("parfile","par");                                 #
# Revision History:                                               #
#         Feb 11 1997    Lianxing Wen Initial Revision            #
#         Feb 23 1997    Lianxing Wen add split(/ /, $c[1])       #
###################################################################
sub findpar{
    my ($parfile, $par) = @_;
    my $par_val;
    if(-f $parfile){
	open(PT,$parfile);
	@l=<PT>;
	foreach(@l){
	    chomp(@l);
	    @c = split("=");
	    if($c[0] eq $par){
		@val = split(/ /,$c[1]);
		$par_val = $val[0];
             }
        }
    }

   close(PT);
    return ($par_val);
}





