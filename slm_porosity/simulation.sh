#!/bin/bash

solver=$(ls wlt_3d_*_heat_transfer.out)
inp=problem.inp
log=log.wlt
res=results
dir=$(pwd)

cases="100 200 400 800 1600 3200"
if test $# -gt 0; then
    cases=$@
fi

for c in $cases; do
    mkdir -p "$dir/$c"
    cd "$dir/$c"
    [ -f $log ] && continue
    echo "Simulate for $c..."
    ( cd ../ && cp $solver $inp $c)
    mkdir $res
    speed=$(grep scanning_speed $inp | awk '{ print $3 }')
    speed=$(python -c "print $c*$speed")
    echo "scanning_speed = $speed" >> $inp
    tail -1 $inp
    ./$solver $inp > $log 2>&1 &
done

