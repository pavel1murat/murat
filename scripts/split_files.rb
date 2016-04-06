#!/usr/bin/env ruby
#-----------------------------------------------------------------------
#  call:     split_dataset.rb fn1 
#
#  assume that fn1 and fn2 have the same structure:
#  4 columns: run_number, than 3 integers
#  merge them into a single file leaving run number only once
#
# used by the lost run sections monitor
#
#  example: 
#  --------
#  cdfopr/scripts/concatenate.rb -d -p 2004_08 -s 11:30
#------------------------------------------------------------------------------
require 'find'
require 'fileutils'
require 'getoptlong'

def usage
  puts "usage: split_dataset --input fn1 --nfiles nfiles_per_fileset"
  exit(-1)
end

usage if ARGV.length < 1

opts = GetoptLong.new(
  [ "--input"   , "-i",        GetoptLong::REQUIRED_ARGUMENT ],
  [ "--nfiles"  , "-n",        GetoptLong::REQUIRED_ARGUMENT ],
  [ "--verbose" , "-v",        GetoptLong::NO_ARGUMENT       ]
)

$dataset = nil;
$nfiles  = 100;
$verbose = 0;

#----------------------------- process the parsed options
opts.each do |opt, arg|
  if    (opt == "--input"         ) ; $dataset = arg
  elsif (opt == "--nfiles"        ) ; $nfiles  = arg.to_i
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
input_file  = File.open($dataset) ; 

fileset     = 1
output_fn   = $dataset+format(".%04i",fileset);
# output_fn   = 'aaaa'+format(".%04i",fileset);
output_file = File.new(output_fn,"w");

nlines = 0;

input_file.each_line { |line| 

  output_file.write(line)
  nlines = nlines + 1
  if (nlines == $nfiles) then
    output_file.close()
    fileset = fileset+1;
#    output_fn   = 'aaaa'+format(".%04i",fileset);
    output_fn   = $dataset+format(".%04i",fileset);
    output_file = File.new(output_fn,"w");
    nlines = 0
  end
}

if (nlines < $nfiles) then
    output_file.close()
end

#-----------------------------------------------------------------------
#  done
#-----------------------------------------------------------------------
puts "split_files.rb: done"
# exit(rc)

