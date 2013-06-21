#!/bin/bash
gcc -c sam2bam.c -pthread -o sam2bam.o
gcc -g -Wall -O2 -o sam2bam sam.o sam2bam.o -Lbcftools libbam.a -lz -lpthread
