#!/bin/bash

  today=`date +%Y-%m-%d`
pattern=$today"-build"
      n=`ls results | grep $today"-build" | wc -l`
   next=$((n+1))
logfile=results/$today-build-`printf "%02i" $next`.log

echo $logfile

scons -j4 >| $logfile 2>&1 &

xterm -g 256x20 -e tail -f $logfile &
