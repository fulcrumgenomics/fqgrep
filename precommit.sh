#!/usr/bin/env bash

###############################################################################
# Script that should be run pre-commit after making any changes.
###############################################################################

set -e

failures=""

function banner() {
    echo
    echo "================================================================================"
    echo $*
    echo "================================================================================"
    echo
}

#####################################################################
# Takes two parameters, a "name" and a "command". 
# Runs the command and prints out whether it succeeded or failed, and
# also tracks a list of failed steps in $failures.
#####################################################################
function run() {
    local name=$1
    local cmd=$2

    banner "Running $name"
    set +e
    $cmd
    exit_code=$?
    set -e
    
    if [[ $exit_code == 0 ]]; then
        echo Passed $name 
    else
        echo Failed $name
        if [ -z "$failures" ]; then
            failures="$name"
        else
            failures="$failures, $name"
        fi
    fi
}

run "cargo fmt"    "cargo fmt --all"
run "cargo clippy" "cargo clippy --all-features --locked -- -D warnings"
run "cargo test"   "cargo test --locked --quiet"

if [ -z "$failures" ]; then
    banner "Precommit Passed"
else
    banner "Precommit Failed with failures in: $failures"
    exit 1
fi
