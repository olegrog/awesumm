#!/bin/bash

server=hilbert

dir=$(pwd);
dir=${dir#$HOME/}

dest=3d/results
[[ $# > 0 ]] && dest=$1

mkdir -p $HOME/$dir/$dest
rsync -auvz --delete $server:$dir/results/*.vtk $HOME/$dir/$dest/
