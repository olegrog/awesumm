#!/bin/bash -e

cases=$(ls | grep '^[0-9]' | sort -n)
if test $# -gt 0; then
    cases=$@
fi

date
printf "%9s %9s\n" 'speed(mm/s)' 'width(um)'
for d in $cases; do
    d=$(basename $d)
    width=$(awk '!/^#/ { print $6 }' $d/results/melting_pool.txt | tail -1)
    width=$(python -c "print $width*5.500e+01")
    printf "%9s %.3e\n" $d $width
done
