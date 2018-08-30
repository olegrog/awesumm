#!/bin/bash

solver=$(ls wlt_3d_*_heat_transfer.out)
res2vtk=$(ls res2vis_*_heat_transfer.out)
inp=problem.inp
log=log.wlt
log2=log.res
res=results
dir=$(pwd)

cases="150 300 600 1200 2400"
if test $# -gt 0; then
    cases=$@
fi

for c in $cases; do
    mkdir -p "$dir/$c"
    cd "$dir/$c"
    [ -f $log ] && continue
    echo "Simulate for $c..."
    ( cd ../ && cp $solver $res2vtk $inp $c)
    mkdir $res
    speed=$(grep scanning_speed $inp | awk '{ print $3 }')
    speed=$(python -c "print $c*$speed")
    length=$(python -c "print 5+int(5*($c/150.)**.5)")
    Mx=$(python -c "print int(5*($c/150)**.25)")
    echo "scanning_speed = $speed" >> $inp
    echo "coord_max = $length, 5.0, 0.0" >> $inp
    echo "M_vector = $Mx, 2, 2" >> $inp
    tail -3 $inp
    ( ./$solver $inp > $log 2>&1; ./$res2vtk results/res.0000.res -t 10000 1 > $log2 2>&1 ) &
done

