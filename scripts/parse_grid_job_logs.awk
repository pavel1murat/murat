#!/usr/bin/awk
#------------------------------------------------------------------------------
#
#  awk -f murat/scripts/parse_g4_log.awk pet_001.log
#------------------------------------------------------------------------------
BEGIN {
    file = fn;
    rc   = -1
    n    = 1;
}

{ 
    if (n == 0)  { 
	printf "  rc  nevents  CPUTime RealTime  file\n";
	n = 1;
    }
}

/TrigReport Events total/ {
    nevents = $5
#   print $5;
}

/TimeReport CPU\/event =/ {
    cpu  = $4;
    real = $7
}

/mu2egrid exit status/ {
    rc = $4
}

END {
    if (rc < 0) {
	nevents = -1
	cpu     = -1
	real    = -1
    }

#    printf "rc: %5i nevents: %6i CPUTime: %12.6f RealTime: %12.6f file: %s", rc, "nevents ", nevents,"CPUTime ", cpu, "RealTime: ", real, fn
    printf "%5i %6i %12.6f %12.6f %s\n", rc,nevents,cpu,real,fn
}
