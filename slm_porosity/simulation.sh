#!/bin/bash

solver=$(ls wlt_3d_*_heat_transfer.out)
res2vtk=$(ls res2vis_*_heat_transfer.out)
inp=problem.inp
res=results
dir=$(pwd)
mpi=1

if [[ $mpi -eq 1 ]]; then
    mpirun="mpirun -np 4"
    com=".com"
fi

cases="150 300 600 1200 2400"
if test $# -gt 0; then
    cases=$@
fi

for c in $cases; do
    mkdir -p "$dir/$c"
    cd "$dir/$c"
    [ -f log.wlt ] && continue
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
    (
        $mpirun ./$solver $inp > log.wlt 2>&1
        grep 'cfl=' log.wlt | nl > log.cfl
        ./$res2vtk results/res.0000$com.res -t 10000 1 > log.res 2>&1
    ) &
done

