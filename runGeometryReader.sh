#!/bin/bash

make

OFFSET=27
for i in {9..12}; do
    # echo "./geometryReader.out FH${i} $((i+OFFSET))"
    ./geometryReader.out FH${i} $((i+OFFSET))
done
    
OFFSET=39
for i in {1..12}; do
    # echo "./geometryReader.out BH${i} $((i+OFFSET))"
    ./geometryReader.out BH${i} $((i+OFFSET))
done
