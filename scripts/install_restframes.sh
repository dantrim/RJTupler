#!/bin/bash

#
# A script to checkout the RestFrames source code and
# compile it if the user wishes to do so.
#
# The checkout and compilation will happen in the caller's
# current working directory.
#
# An analysis base should be set up (primarily for obtaining
# the ROOT version that will be associated with the analysis
# code that is using RestFrames.
#
# daniel.joseph.antrim@cern.ch
# October 2018
#

analysis_base="21.2.45"
default_tag="master"

function print_usage {
    echo ""
    echo "This script will checkout the RestFrames source code and"
    echo "compile it in the user's current working directory."
    echo ""
    echo "usage: source <path>/install_restframes.sh [options]"
    echo "options:"
    echo " --no-compile             : do not perform the compilation [default: false]"
    echo " -t|--tag                 : choose a specific tag of RestFfames [default: $default_tag]"
    echo " -h|--help                : print this help message"
    echo ""

}

function install_rf {

    tag=${1}
    no_compile=${2}

    asetup "AnalysisBase,${analysis_base}"

    git clone https://github.com/crogan/RestFrames RestFrames
    if [[ ! -d "RestFrames" ]]; then
        echo "install_restframes    ERROR Failed to checkout RestFrames from git!"
        return 1
    fi

    pushd RestFrames/
    ./configure --prefix=${PWD}
    make -j4 2>&1 |tee compile.log
    make install -j4 2>&1 |tee install.log
    popd

}

function main {

    tag=${default_tag}
    no_compile="0"

    while test $# -gt 0
    do
        case $1 in
            --no-compile)
                no_compile"1"
                ;;
            -t)
                tag=${2}
                shift
                ;;
            --tag)
                tag=${2}
                shift
                ;;
            -h)
                print_usage
                return 0
                ;;
            --help)
                print_usage
                return 0
                ;;
            *)
                echo "install_restframes    ERROR Invalid argument: ${1}"
                return 1
                ;;
        esac
        shift
    done

    echo "-----------------------------------------------"
    echo "install_restframes   start          : `date`"
    echo "install_restframes   no compilation : ${no_compile}"
    echo "install_restframes   restframes tag : ${tag}"
    echo "-----------------------------------------------"

    if ! install_rf $tag $no_compile ; then
        echo "install_restframes    failed"
        return 1
    fi
    echo "install_restframes    success"
    


}

#________
main $*
