#
#  awk -f murat/scripts/parse_g4_log.awk pet_001.log
#------------------------------------------------------------------------------
BEGIN{
    n=0
}

{ 
    if (n == 1) print $0; 
    n = 0;
}

/Step#/ {
    if ($1 != "#Step#") print $0;
    n = 1;
}

/G4Track Information:/ {
    print "-------------------------------------------------------------------------------------------------------"
    print $0
}

END {
}
