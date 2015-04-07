#!/usr/bin/env ruby
#-----------------------------------------------------------------------
#  call:     find_missing_events.rb fn1 fn2 
#
#  assume that fn1 and fn2 have the same structure:
#  1 line per event
#-----------------------------------------------------------------------
require 'find'
require 'fileutils'

def usage
  puts "usage: find_missing_events fn1 fn2"
  exit(-1)
end

usage if ARGV.length < 2

fn1     = ARGV[0]
fn2     = ARGV[1]
rc          = 0;
#-----------------------------------------------------------------------
#  start 
#-----------------------------------------------------------------------
f1          = File.open(fn1) ; 

f1.each_line { |line1| 

   ok = 0

   f2 = File.open(fn2) ; 
   f2.each_line { |line2| 

     if (line1 == line2) 
       ok = 1
       break
     end
   }

   if (ok == 0) 
       printf "missing in %s : %s", fn2, line1
   end
   f2.close;
}
#-----------------------------------------------------------------------
#  done
#-----------------------------------------------------------------------
f1.close;

exit(rc)

