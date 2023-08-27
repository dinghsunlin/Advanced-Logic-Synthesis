#!/bin/bash

DIR=$PWD
BENCHMARK="./benchmarks"
LOG="/users/student/mr111/dhlin22/ALS/Final/log"
RESULT="/users/student/mr111/dhlin22/ALS/Final/results"

SISBIN="/users/student/mr111/dhlin22/ALS/Final/SIS/bin"

EXE="/users/student/mr111/dhlin22/ALS/Final/Power_Optimization_State_Assignment"
TMP="/users/student/mr111/dhlin22/ALS/Final/tmp"
SIS="/users/student/mr111/dhlin22/ALS/Final/SIS/bin/sis"
SOURCE="/users/student/mr111/dhlin22/ALS/Final/SIS/bin/opt_map_power.scr"
POWER="/users/student/mr111/dhlin22/ALS/Final/Power_Results.log"

make clean
make

echo -n "" > $POWER

for file in $BENCHMARK/*
do
    bench=$(echo $file | cut -d "/" -f 3)
    bench=$(echo $bench | cut -d "." -f 1)
    echo $bench
    $EXE $file > $LOG/$bench.log
    echo "" >> $LOG/$bench.log

    line="read_blif "$RESULT"/"$bench".blif"
    echo $line > $TMP
    line="source "$SOURCE
    echo $line >> $TMP
    echo "quit" >> $TMP
    
    cd $SISBIN
    $SIS -f $TMP -x >> $LOG/$bench.log
    $SIS -f $TMP -x >> $POWER
    echo "" >> $POWER
    cd $DIR
    rm -rf $TMP
done
