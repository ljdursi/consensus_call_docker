#!/bin/bash

readonly COMMAND=$1
readonly SUBCOMMAND=$2

readonly EXECUTABLE_PATH=${USE_EXECUTABLE_PATH:-"/usr/local/bin"}

function usage {
    echo >&2 "$0: execute a command from the docker. "
    echo >&2 "    options: "
    echo >&2 "       consensus snv [opts]: make consensus SNV calls"
    echo >&2 "       consensus indel [opts]: make consensus indel calls"
    echo >&2 "       filter {type} [opts]: filter consensus calls"
    echo >&2 "       download reference /dbs/path: download, install pancan reference"
    echo >&2 "       download annotations /dbs/path: download, install annotations"
    exit 1
}

###
### Make sure options are valid 
###
if [[ -z "${COMMAND}" ]] || [[ -z "${SUBCOMMAND}" ]] 
then
    echo >&2 "Missing arguments."
    usage
fi

if [[ "${COMMAND}" != "consensus" ]] && [[ "${COMMAND}" != "filter" ]] && [[ "${COMMAND}" != "download" ]]
then
    echo >&2 "Invalid command ${COMMAND}."
    usage
fi

###
### if $COMMAND is consensus, run consensus step:
###
if [[ "${COMMAND}" == "consensus" ]] 
then
    if [[ "${SUBCOMMAND}" != "snv" ]] && [[ "${SUBCOMMAND}" != "indel" ]] 
    then
        echo >&2 "Invalid variant type ${SUBCOMMAND}"
        usage
    fi

    "${EXECUTABLE_PATH}/consensus_${SUBCOMMAND}.sh" "${@:3}"
fi

###
### if $COMMAND is download, download the data
###
if [[ "${COMMAND}" == "download" ]] 
then
    "${EXECUTABLE_PATH}/build_dbs.sh" "${@:2}"
fi

###
### if $COMMAND is filter, run the appropriate filter step
###
if [[ "${COMMAND}" == "filter" ]] 
then
    "${EXECUTABLE_PATH}/filter/additional_filters.sh" "${@:2}"    
fi
