#!/bin/sh
#Helper script to filter and analyze a given category
Category=$1

if [ -z "${Category}" ]; then
    echo "Usage: $0 Category CategoryValues"
    exit 1
fi