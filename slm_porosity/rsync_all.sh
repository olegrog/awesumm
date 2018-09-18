#!/bin/bash

srv=hil

dir=$(pwd)
dir=${dir#$HOME/}
dest=3d
[[ $# > 0 ]] && dest=$1

mkdir -p $HOME/$dir/$dest
cases=$(ssh $srv "cd $dir; ls | grep ^[0-9]")

for c in $cases; do
    echo "Copying $c to $dest..."
    rsync -auvz $srv:$dir/$c/results/*.{vtk,txt} $HOME/$dir/$dest/$c/
    rsync -auvz $srv:$dir/$c/*.inp $HOME/$dir/$dest/$c/
    rsync -auvz $srv:$dir/$c/log.cfl $HOME/$dir/$dest/$c/
done
