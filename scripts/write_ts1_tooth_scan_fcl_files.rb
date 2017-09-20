#!/usr/bin/env ruby
#-----------------------------------------------------------------------
#  call:     write_ts1_tooth_scan_fcl_files.rb -o ../grid_fcl/ts1_tooth
#
#  example: 
#  --------
#  cdfopr/scripts/write_scan.rb 
#------------------------------------------------------------------------------
require 'find'
require 'fileutils'
require 'getoptlong'

def usage
  puts "usage: write_scan --output dir"
  exit(-1)
end

# usage if ARGV.length < 1

opts = GetoptLong.new(
  [ "--output"  , "-o",        GetoptLong::REQUIRED_ARGUMENT ],
  [ "--nevents" , "-n",        GetoptLong::REQUIRED_ARGUMENT ],
  [ "--verbose" , "-v",        GetoptLong::NO_ARGUMENT       ]
)
$nevents = 1
$odir    = "junk";
$verbose = 0;
#----------------------------- process the parsed options
opts.each do |opt, arg|
  if    (opt == "--output"        ) ; $odir    = arg
  elsif (opt == "--nevents"       ) ; $nevents = arg;
  elsif (opt == "--verbose"       ) ; $verbose = 1;
  end
 
  if ($verbose != 0) ; puts "Option: #{opt}, arg #{arg.inspect}" ; end
end

# pwd=`pwd`.strip
# puts "emoe input TCL file = #{$input_tcl_file}"
# puts "smin = #{$smin} smax = #{$smax}"
#-----------------------------------------------------------------------
#  start 
#-----------------------------------------------------------------------

$height = [ 5,10,20,30,40,50] ;
$width  = [50,60,70,80,90,100, 150, 200] ;
$length = [20,30,40,50,80,100,150,200];

for ih in 0..$height.length-1
  height = $height[ih];
#------------------------------------------------------------------------------
# different subdirectories for different heights 
# for each height the job numbering starts from zero
# reserve "622_0000" for the "no tooth" configuration
#------------------------------------------------------------------------------
  job = 0;
  od1 = $odir+"/"+format("622_%04i",ih+1);   # $dataset+format(".%04i",fileset);

  for iw in 0..$width.length-1
    width  = $width[iw];
    for il in 0..$length.length-1
      length = $length[il];

      output_dir = od1+"/"+format("%03i_%03i_%03i",job,iw,il);   # $dataset+format(".%04i",fileset);
      `mkdir -p #{output_dir}`
      `cp murat/test/ts1_tooth/622_0000/* #{output_dir}`
#------------------------------------------------------------------------------
# modify the geometry definition file
#------------------------------------------------------------------------------
      fn_geom=output_dir+"/geom.txt"
      f_geom = File.open(fn_geom) 
      fn_geom_2 = fn_geom+"2";
      f_geom_2 = File.new(fn_geom_2,"w");
      f_geom.each_line { |line| 
#        puts line
        if (line.index("ts.coll1.useFlashBlock") == nil) then
          f_geom_2.write(line)
        end
      }
      f_geom.close();

      f_geom_2.write("bool   ts.coll1.useFlashBlock      =  true ;\n");
      f_geom_2.write(format("double ts.coll1.flashBlock.Height  = %5.1f ; // mm ;\n",height.to_f));
      f_geom_2.write(format("double ts.coll1.flashBlock.Width   = %5.1f ; // mm ;\n",width.to_f));
      f_geom_2.write(format("double ts.coll1.flashBlock.Length  = %5.1f ; // mm ;\n",length.to_f));
      f_geom_2.close();
      `mv #{fn_geom_2} #{fn_geom}`
#------------------------------------------------------------------------------
# update beam_g4s1.fcl
#------------------------------------------------------------------------------
      fn1 = output_dir+"/beam_g4s1.fcl"
      f1  = File.open(fn1) 
      fn2 = fn1+".tmp";
      f2  = File.new(fn2,"w");
      f1.each_line { |line| 
#        puts line
        if (line.index("services.GeometryService.inputFile") == nil) then
          f2.write(line)
        else
          new_line=format("services.GeometryService.inputFile              : \"%s/geom.txt\"\n",output_dir);
          f2.write(new_line);
        end
      }
      f2.write("source.firstRun            : 1\n");
      f2.write("source.firstSubRun         : 1\n");
      f2.write(format("source.maxEvents           : %i\n",$nevents));
      f2.write("mu2emetadata.fcl.prologkeys: [  ]\n");
      f2.write("mu2emetadata.fcl.inkeys    : [  ]\n");
      f2.write("mu2emetadata.fcl.outkeys   : [  ]\n");

      f1.close();
      f2.close();
      `mv #{fn2} #{fn1}`
#------------------------------------------------------------------------------
# finally, increment the job number 
#------------------------------------------------------------------------------
      job = job +1;
    end
  end
end
#-----------------------------------------------------------------------
#  done
#-----------------------------------------------------------------------
puts "write_scan.rb: done"
# exit(rc)

