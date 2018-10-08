#!/bin/bash


function main {

    cwd=${PWD}
    ana_dir="susynt-read"
    if [[ ! ${cwd} == *"${ana_dir}"* ]]; then
        echo "setup_restframes    ERROR susynt-read/ is not in current working directory!"
        return 1
    fi

    path_to_read=$(echo $cwd | sed 's/\(susynt-read\).*/\1/g')

    rj_dir=${path_to_read}/RestFrames
    if [[ ! -d $rj_dir ]]; then
        echo "setup_restframes    ERROR RestFrames/ directory is not found in expected location (=${rj_dir})"
        return 1
    fi

    source ${rj_dir}/setup_RestFrames.sh
    export CMAKE_PREFIX_PATH=${rj_dir}/lib:${CMAKE_PREFIX_PATH}

}

main $*
