#!/bin/bash
# git clone https://github.com/soedinglab/hh-suite
export REFORMAT_SCRIPT_PATH="$HOME/hh-suite/scripts"
# $1: folder containing .a3m files to expand
ali_path="$1"

# convert a3ms to a2ms if needed, including insertions in a2m output
# skip a2ms that are already present and have the correct number of sequences
reformat_a3m() {
    a3mfile="$1"
    a2mfile="${a3mfile/.a3m/.a2m}"
    if [[ -f $a2mfile && -s $a2mfile ]] && (( "$(grep -c '>' $a3mfile)" == "$(grep -c '>' $a2mfile)" )); then
        return
    fi
    echo "reformat.pl $a3mfile -> $a2mfile"
    $REFORMAT_SCRIPT_PATH/reformat.pl a3m a2m $a3mfile $a2mfile
}
export -f reformat_a3m
# run num-cpus number of jobs
ls "$ali_path"/*.a3m | xargs -n 1 -P "$(nproc)" -I {} bash -c 'reformat_a3m "$@"' _ {}
