#!/bin/bash
# wget http://eddylab.org/software/hmmer/hmmer.tar.gz && tar zxf hmmer.tar.gz && cd hmmer-3.3.2 && ./configure --prefix ~/hmmer && make -j$(nproc) && make install
export HMMER_PATH="$HOME/hmmer"
ali_path="$1"

# build HMMs
build_hmm() {
    alifile="$1"
    hmmfile="${alifile/.a2m/.hmm}"
    if [[ -s $hmmfile ]]; then
        continue
    fi
    echo "hmmbuild $hmmfile $alifile"
    $HMMER_PATH/bin/hmmbuild --amino --hand "$hmmfile" "$alifile"
}
export -f build_hmm
ls "$ali_path"/*.a2m | xargs -n 1 -P "$(nproc)" -I {} bash -c 'build_hmm "$@"' _ {}
