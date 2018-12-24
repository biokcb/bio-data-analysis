#!/bin/bash

# Download the sequencing data using SRA toolkit
for i in $( cat samples.txt ); do
    fasterq-dump $i
done