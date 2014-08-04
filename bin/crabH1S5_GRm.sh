#!/bin/bash

NINST=$1
NINJ=$2

arg=("J0534+2200" "H1" "S5" "GR" "m" "$NINST" "$NINJ")

./injsrch_master.py "${arg[@]}"

condor_submit subs/injsrch_H1S5_J0534+2200_GRm.sub

