#!/bin/bash

PSR=$1
KIND=$2
NINST=$3
NINJ=$4

arg=("$PSR" "H1" "S5" "$KIND" "p" "$NINST" "$NINJ")

./injsrch_master.py "${arg[@]}"

condor_submit subs/injsrch_H1S5_"$PSR"_"$KIND"p.sub

