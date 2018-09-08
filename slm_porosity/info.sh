#!/bin/bash -e

cases=$(ls | grep '^[0-9]' | sort -n)
if test $# -gt 0; then
    cases=$@
fi

date
printf "%9s %9s %9s\n" 'speed(mm/s)' 'width(um)' 'progress(%)'
for d in $cases; do
    d=$(basename $d)
    width=$(awk '!/^#/ { print $6 }' $d/results/melting_pool.txt | tail -1)
    width=$(python -c "print $width*2.750e+01")
    if test -f $d/log.wlt; then
        time=$(awk '/^t= / { print $2 }' $d/log.wlt | tail -1)
        len=$(awk '/coord_max = / { print $3 }' $d/problem.inp | sed 's/,//' | tail -1)
        speed=$(awk '/scanning_speed = / { print $3 }' $d/problem.inp | tail -1)
        progress=$(python -c "print 100*$time/($len-5.)*$speed")
    fi
    printf "%9s %.3e %.1f\n" $d $width $progress
done
