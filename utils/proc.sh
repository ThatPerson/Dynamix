#!/bin/bash

cat $1/orderparams_proc* > $1/orderparams.dat
cat $1/final_proc* > $1/final.dat
cat $1/errors_proc* > $1/errors.dat
cat $1/gaf_proc* > $1/gaf.dat
