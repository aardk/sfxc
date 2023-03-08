#!/usr/bin/perl 
#===========================================================================
#
# Copied from DiFX "geteop.pl 4721 2012-07-18 15:55:03Z HelgeRottmann"
#
# Adapted by Jan Wagner for SFXC. Compared to DIFX, SFXC expects a 
# different format for the VEX $EOP; section -- a single multi-datapoint 
# EOP def. instead of DiFX's multiple single-datapoint EOP def's.
#
#============================================================================

$FINALS_FILE = "usno_finals.erp";
#$FINALS_URL = "http://gemini.gsfc.nasa.gov/solve_save/$FINALS_FILE";
#$FINALS_URL = "ftp://ftp.lbo.us/pub/staff/wbrisken/EOP/$FINALS_FILE";
#$FINALS_URL = "ftp://cddis.gsfc.nasa.gov/vlbi/gsfc/ancillary/solve_apriori/$FINALS_FILE";
$FINALS_URL = "ftp://gdc.cddis.eosdis.nasa.gov/vlbi/gsfc/ancillary/solve_apriori/$FINALS_FILE";
$EOPVEX_FILE = "EOP_SFXC.txt";

sub printUsage
{
        print "----------------------------------------------------------------------------------\n";
        print " Script to download and extract EOP solutions from URL\n";
	print "   $FINALS_URL\n";
	print " The EOPS are reformatted so they can be used in a VEX file for SFXC correlation.\n";
        print "---------------------------------------------- -----------------------------------\n";
        print "\nUsage: $0 yyyy-doy numEOP\n\n";
        print "  yyyy-doy: The year and day-of-year of the first daily EOP value to obtain\n";
        print "            This date may be chosen to be 2-3 days before the VLBI experiment.\n"; 
        print "  numEOP:   The number of EOPs to get. Note that SFXC has a limit of 10 EOPs.\n\n";
        print "The output will be written to $EOPVEX_FILE.\n\n";
        print "---------------------------------------------- -----------------------------------\n";
	exit;
}


$eopCount = 0;

if ($#ARGV != 1)
{
	&printUsage;
}
	
$date = $ARGV[0];
$numEop = $ARGV[1];

# split date into year and doy, calculate julian date
if ($date =~ /(\d{4})-(\d+)/)
{
	$year = $1;
	$doy = $2;
	#print "year: $year $doy\n";
        #$jdate = &calcJD($year,$doy, 23, 59, 59);
        $jdate = &calcJD($year,$doy, 0, 0, 0);
}
else
{
	&printUsage;
}

if ($jdate >= 2457754.5) 
{
	# Leap second 01.Jan 2017
	$TAIUTC = 37;
}
elsif ($jdate > 2457203.5) 
{
	# Leap second 30.Jun 2015
	$TAIUTC = 36
}
elsif ($jdate > 2456109.5) 
{
	# Leap second 30.Jun 2012
	$TAIUTC = 35;
}
else
{
	$TAIUTC = 34;
}

# get a recent USNO solution file
if ( -f $FINALS_FILE )
{
	my $now = time();
	my @fstats = stat($FINALS_FILE);
	if ($now-$fstats[9] > 2*86400) {
		system ("rm $FINALS_FILE");
                system ("curl -u anonymous:jops\@jive.nl --ftp-ssl $FINALS_URL -o $FINALS_FILE");
	} else {
		print "Using existing $FINALS_FILE (less than 2 days old).\n";
	}
}
else
{
        system ("curl -u anonymous:jops\@jive.nl --ftp-ssl $FINALS_URL -o $FINALS_FILE");
}
sleep(1);


open(INFILE, $FINALS_FILE) || die("Could not open file: $FINALS_FILE");
open(OUTFILE, ">$EOPVEX_FILE") || die("Could not open file: $EOPVEX_FILE");

print OUTFILE  "\$EOP;\n";
print OUTFILE  "  def EOP_SFXC;\n";
print OUTFILE  "    TAI-UTC        = $TAIUTC sec;\n";
print OUTFILE  "    A1-TAI         = 0 sec;\n";
print OUTFILE  "    num_eop_points = $numEop;\n";
print OUTFILE  "    eop_interval   = 24 hr;\n";

$eopDay0 = 0;
my $ut1utcList = "";
my $xwobList = "";
my $ywobList = "";

# parse the file
while (<INFILE>)
{
	# skip comment lines
	if ($_ =~ /#/)
	{
		next;
	}

	@fields = split(/\s+/, $_);

	# check if first column contains a valid JD
	if ($fields[0] =~ /\d+\.\d+/)
	{	

		if (($fields[0] ge $jd) && ($eopCount < $numEop))
		{
			#print "$eopCount $numEop\n";
			#print "0:$fields[0] 1:$fields[1] 2:$fields[2] 3:$fields[3]\n";
                        print "t = $fields[0], jd = $jd \n";
			$xWobble = $fields[1] / 10.0;
			$yWobble = $fields[2] / 10.0;
			$ut1utc  = $fields[3] / 1e6 + $TAIUTC; 
			$eopDay = $doy + $eopCount;
			if ($eopDay0 == 0)
			{
				$eopDay0 = $eopDay
			}

			$ut1utcList = $ut1utcList . sprintf(" %.6f sec :", $ut1utc);
			$xwobList = $xwobList . sprintf(" %.6f asec :", $xWobble);
			$ywobList = $ywobList . sprintf(" %.6f asec :", $yWobble);

			$eopCount += 1;
		}
	}
}

# tidy up, remove last ':'
chop($ut1utcList);
chop($xwobList);
chop($ywobList);

# create output
print OUTFILE "    eop_ref_epoch=$year" . "y$eopDay0" . "d;\n";
printf OUTFILE ("    ut1-utc  = %s ;\n", $ut1utcList);
printf OUTFILE ("    x_wobble = %s ;\n", $xwobList);
printf OUTFILE ("    y_wobble = %s ;\n", $ywobList);
print OUTFILE "  enddef;\n";

print "Output written to $EOPVEX_FILE: \n";
system ("cat $EOPVEX_FILE");

exit;

sub calcJD
{
        my $year = $_[0];
        my $day = $_[1];
        my $hour = $_[2];
        my $minute = $_[3];
        my $second = $_[4];

        $y1= $year - 2000;
        $im = 10;
        $iy = int($year -1);
        $ic = int($iy /100);
        $iy = int($iy - $ic *100);
        $ix1 = int(146097 * $ic /4);
        $ix2 = int(1461 * $iy /4);
        $ix3 = int((153 * $im + 2) /5);
        $jd= $ix1 + $ix2 + $ix3 + $day - 678882;

        $jd += $hour/24 + $minute/1440 + $second / 86400 + 2400000.5;

        return($jd)
}
