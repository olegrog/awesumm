#!/bin/bash -e

cases=$(ls | grep '^[0-9]' | sort -n)
if test $# -gt 0; then
    cases=$@
fi

time_stat() {
    if [ "$(uname)" == "Darwin" ]; then
        stat -l -t '%s' $1 | awk '{print $6}' | sort
    else
        stat --format='%Z' $1 | sort
    fi
}

humanize() {
    local T=$1
    local D=$((T/60/60/24))
    local H=$((T/60/60%24))
    local M=$((T/60%60))
    local S=$((T%60))
    [[ ! $T ]] && return
    [[ $D > 0 ]] && printf ' %dd' $D
    [[ $H > 0 ]] && printf ' %2dh' $H
    [[ $M > 0 ]] && printf ' %2dm' $M
    [[ $D = 0 && $H = 0 ]] && printf ' %ds' $S
}

fmt=$(printf '%%15s %.0s' {1..5})
date
printf "$fmt\n" 'speed(mm/s)' 'width(um)' 'progress(%)' 'spend_time' 'remaining_time'
for d in $cases; do
    d=$(basename $d)
    width=$(awk '!/^#/ { print $6 }' $d/results/melting_pool.txt | tail -1)
    [[ -z $width ]] && width=0
    width=$(python -c "print($width*2.750e+01)")
    if test -f $d/log.wlt; then
        time=$(awk '/^t= / { print $2 }' $d/log.wlt | tail -1)
        len=$(awk '/coord_max = / { print $3 }' $d/problem.inp | sed 's/,//' | tail -1)
        speed=$(awk '/scanning_speed = / { print $3 }' $d/problem.inp | tail -1)
        progress=$(python -c "print(100*$time/($len-5.)*$speed)")
    fi
    first_time=$(time_stat $d/problem.inp)
    last_time=$(time_stat $d/log.wlt)
    spend_time=$(($last_time-$first_time))
    remaining_time=$(python -c "print(int($spend_time*(100-$progress)/$progress))")
    width=$(printf %.3e $width)
    progress=$(printf %.1f $progress)
    printf "$fmt\n" $d $width $progress "$(humanize $spend_time)" "$(humanize $remaining_time)"
done
