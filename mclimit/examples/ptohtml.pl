#!/usr/bin/perl

$found=0;
while($line=<>)
{
    if ($line =~ /csm_model::print  -- printing out model information/) { $found=1; last; }
}
exit if ($found == 0);

print << 'EOF';
<HTML>
<HEAD>
<TITLE>Channel Systematic Summary</TITLE>
</HEAD>
EOF

while($line=<>)
{
    #print "Input line: ",$line;
    last if ( ($line =~ /-------------------/) && $numchans > 0);
    if ($line =~ /Channel:/) {
	chomp($line);
	$line =~ /Name: /;
	$channame=$';
	print "<h3>Channel name: ",$channame,"</h3>\n";
        procchannel();
    }
}

print << 'EOF';
</HTML>
EOF

exit;


sub procchannel {

    my @sys = ();
    my @histoname = ();
    my @histotitle = ();
    my %snhash = ();
    my @integral = ();
    my @sigflag = ();
    my @poissflag = ();
    my $sigtot=0;
    my $bgtot=0;

    while ($line2=<>) {
	#print "line2 input: ",$line2;
        if ($line2 =~ /End/)
	{ #print "end of channel input.\n";
	  last;}
	chomp($line2);
	if ($line2 =~ /bbeta/) {}
	elsif ($line2 =~ /Template /) {
	    $templatenum = $';
	    #print "Found a template: <",$templatenum,">\n";
	}
	if ($line2 =~ /Histogram name:/) {
	    $l3p1 = $';
	    $l3p1 =~ /Histogram title:/;
	    $histoname[$templatenum] = $`;
	    $l3p2 = $';
	    $l3p2 =~ /sft:/;
	    $histotitle[$templatenum] = $`;
	}
	if ($line2 =~ /signalflag: /) {
	    $sigflag[$templatenum] = $';
	}
	if ($line2 =~ /poissflag: /) {
	    $poissflag[$templatenum] = $';
	}
	if ($line2 =~ /Scaled Integral with all syst:/) {
	    $integral[$templatenum] = $';
	    if ($sigflag[$templatenum] =~ /1/) {
		$sigtot += $integral[$templatenum]; 
	    }
	    else {
		$bgtot += $integral[$templatenum]; 
	    }
	}
	if ($line2 =~ /Syst:/) {
	    $sysname = $';
	    $snhash{$sysname} = 1;
	}
	if ($line2 =~ /Up rate error:/) {
	    $sys[$templatenum]{$sysname}{"up"} = $';
	}
	if ($line2 =~ /Down rate error:/) {
	    $sys[$templatenum]{$sysname}{"down"} = $';
	}
	if ($line2 =~ /Total up rate error incl. shape:/) {
	    $sys[$templatenum]{$sysname}{"up"} = $';
	    $sys[$templatenum]{$sysname}{"upshape"} = "s";
	}
	if ($line2 =~ /Total down rate error incl. shape:/) {
	    $sys[$templatenum]{$sysname}{"down"} = $';
	    $sys[$templatenum]{$sysname}{"downshape"} = "s";
	}
    }

    print "<TABLE BORDERCOLOR=BLACK BORDER=2 RULES=ALL FRAME=BOX>\n";
    print "<TR><TD>&nbsp;</TD>";
    for ($i=0; $i < @histotitle; $i++) {
	print "<TD> ",$histoname[$i]," </TD>";
    }
    print "</TR>\n";
    print "<TR><TD>&nbsp;</TD>";
    for ($i=0; $i < @histotitle; $i++) {
	print "<TD> ",$histotitle[$i]," </TD>";
    }
    print "</TR>\n";

    print "<TR><TD>Rate</TD>";
    for ($i=0; $i < @histotitle; $i++) {
	print "<TD> ",$integral[$i]," </TD>";
    }
    print "</TR>\n";

    print "<TR><TD>Signal Flag</TD>";
    for ($i=0; $i < @histotitle; $i++) {
	print "<TD> ",$sigflag[$i]," </TD>";
    }
    print "</TR>\n";

    print "<TR><TD>Poisson Flag</TD>";
    for ($i=0; $i < @histotitle; $i++) {
	print "<TD> ",$poissflag[$i]," </TD>";
    }
    print "</TR>\n";

    foreach $k1 (sort keys %snhash) {
	print "<TR><TD>",$k1," (up)</TD>";
	for ($i=0; $i < @histotitle; $i++) {
	    print "<TD> ",$sys[$i]{$k1}{"up"},$sys[$i]{$k1}{"upshape"}," </TD>";
        }
	print "</TR>\n<TR><TD>",$k1," (down)</TD>";
	for ($i=0; $i < @histotitle; $i++) {
	    print "<TD> ",$sys[$i]{$k1}{"down"},$sys[$i]{$k1}{"downshape"}," </TD>";
        }
	print "</TR>\n";
    }
    print "</TABLE>\n";
    print "<BR>Total Signal: ",$sigtot,"<BR>Total Background: ",$bgtot,"<P><P>\n";

# example working printout

#    for ($i=0; $i < @histotitle; $i++) {
#	print $histoname[$i],"\n";
#	print $histotitle[$i],"\n";
#	foreach $k1 (sort keys %snhash) {
#	    print $k1," up: ",$sys[$i]{$k1}{"up"},$sys[$i]{$k1}{"upshape"},
#	    " down: ",$sys[$i]{$k1}{"down"},$sys[$i]{$k1}{"downshape"},"\n";
#	}
#    }

}
