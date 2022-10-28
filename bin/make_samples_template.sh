#!/bin/bash

set -e
set -u
set -o pipefail

#find -L ../raw_data/ -name '*R1*' | perl -npe 's/\.\.\/raw_data\///' | sort | perl -ne 'BEGIN{print "sample\tcontrol\tfq1\tfq2\tsample_group\n"}; chomp; /([^_\/]+)_[^\/\s]+$/; $samp=$1; print "$samp\tNA\t"; print "$_\t"; s/_R1_/_R2_/; print "$_\t$samp\n"' > samples_template.tsv

module load bbc/R/R-4.1.0-setR_LIBS_USER
Rscript make_samples_template.R 
