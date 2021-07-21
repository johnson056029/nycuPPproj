#!/bin/sh

if [ "$1" = "mine" ]
then
        cp src/viterbi_mine.cpp src/viterbi.cpp
        cp src/trainingSeq_orig.cpp src/trainingSeq.cpp
        cp src/StochHMM_mine.cpp src/StochHMM.cpp
elif [ "$1" = "orig" ]
then
        cp src/viterbi_original.cpp src/viterbi.cpp
        cp src/trainingSeq_orig.cpp src/trainingSeq.cpp 
        cp src/StochHMM_orig.cpp src/StochHMM.cpp
else
        echo "error input!\n"
fi
